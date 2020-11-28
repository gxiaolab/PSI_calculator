#!/usr/bin/python
import sys
import os
import pysam 
from glob import glob 
import glob 
from collections import *
from optparse import OptionParser
from os import getpid, makedirs, path, system # Kind of essential
from time import strftime, time # Not essential
from HTSeq import GenomicInterval as gi
import HTSeq 
import pandas as pd 
import itertools


__author__  = "giovas"
__last_updated__ = "11/27/2020"


"""
---- Update history ---
11/26/2020 Better handling of overlapping exons in annotation
03/08/2019 Fixing of Last exon PSI 
01/08/2019 Correction for NER calculation
11/02/2018 Non strand-specific option bug 
10/25/2018 Include read filters (multi-mapped)
10/12/2018 Replace bedtools with pysam
-----------------------

For non-strand specific data, reads are assumed 
to map to the strand where the annotated exon is

For single-end data  
"""

def print_memory_usage(ExtraMessage = ""):
    message = strftime("[%a, %d %b %Y %H:%M:%S]\t{}".format(ExtraMessage))
    sys.stderr.write("{}\n".format(message)) 

def strands(read, strandedness):
    if strandedness == 'fr-firststrand':
        strand1 = '-' 
        strand2 = '+'
    elif strandedness == 'fr-secondstrand':
        strand1 = '+' 
        strand2 = '-'
    elif strandedness == 'fr-unstrand':
        strand1 = '.'
        strand2 = '.'
    strand = ''
    ## For paired-end (PE)
    if read.is_read1 and not read.is_reverse:
        strand = strand1
    elif read.is_read1 and read.is_reverse:
        strand = strand2
    elif read.is_read2 and read.is_reverse:
        strand = strand1
    elif read.is_read2 and not read.is_reverse: 
        strand = strand2
    ## For single-end (SE)
    else:
        if not read.is_reverse:
            strand = strand1
        elif read.is_reverse:
            strand = strand2
    return strand


def process_CIGAR(cigar):
    match_lengths = [] # Matching portion of the reads
    block_lengths = [] # Genomic coordinates of the Matching portions
    junct_lengths = [] # Length of junctions
    matches      = 0 
    blocks       = 0 
    no_of_blocks = 1
    for (ctype, length) in cigar:
        if ctype == 0: #M
            matches += length
            blocks  += length
        elif ctype == 3: #N
            junct_lengths.append(length)
            match_lengths.append(matches)
            block_lengths.append(blocks)
            blocks += length
            matches = 0 
            no_of_blocks += 1
        elif ctype == 2: #D
            blocks += length
    block_lengths.append(blocks)
    match_lengths.append(matches)
    assert len(block_lengths) == no_of_blocks
    assert len(match_lengths) == no_of_blocks
    return match_lengths, block_lengths, junct_lengths



def read_gtf_file(gtf_file, strandedness):
    exons = list()
    with open(gtf_file, 'r') as gtf_h:
        for line in gtf_h:
            if line.startswith("#"):
                continue
            L = line.strip('\n').split('\t') 
            chrom, source, feature_type, feature_start, feature_end, dot1, strand, dot2, attributes = L
            feature_start, feature_end = map(int, [feature_start, feature_end])

            if feature_type in ['gene', 'transcript']: 
                continue

            original_strand = strand
            
            if strandedness == 'fr-unstrand':
               strand = '.'
            feature_id = attributes.replace(' ', '_')

            feature_length = feature_end - feature_start + 1

            exons.append((chrom, strand, feature_start - 1, feature_end, feature_id, feature_type, feature_length, original_strand, 0, 0, 0.0))

    df = pd.DataFrame(exons, columns = ['CHROM', 'STRAND', "START", "END", "ID", "TYPE", "LENGTH", "STRAND_ORIGINAL", "INCLUSION", "EXCLUSION", "PSI"])

    df = df.sort_values(by = ['CHROM', 'STRAND', 'START', 'END'], ascending = [True, True, True, True])

    df = df.reset_index(drop = True)

    del exons

    return df

      


def get_coverage(exons, bam_file,  overhang):

    bam = pysam.Samfile(bam_file, 'rb')
    read_lengths = defaultdict(int)

    all_chroms  = exons.CHROM.unique()
    all_strands = exons.STRAND.unique()

    for (chrom, fstrand) in itertools.product(all_chroms, all_strands):

        baseline_index = exons[(exons["CHROM"] == chrom) & (exons["STRAND"] == fstrand)].index[0]

        if chrom not in bam.references:
            continue      

        def get_exon(overlap_index):
            if overlap_index >= exons.shape[0]:
                return gi(chrom, sys.maxint - 2, sys.maxint)
            elif exons.iloc[overlap_index].loc["CHROM"] != chrom or exons.iloc[overlap_index].loc["STRAND"] != fstrand:
                return gi(chrom, sys.maxint - 2, sys.maxint)
            else:
                me_start = exons.iloc[overlap_index].loc["START"]
                me_end   = exons.iloc[overlap_index].loc["END"]
                return gi(chrom, me_start, me_end)

        r_count = 0
        ov_index = baseline_index
        overlapping_exon = get_exon(ov_index)

        for read in bam.fetch(chrom, until_eof = True): 
            
            strand = strands(read, strandedness)

            if strand != fstrand:
                continue

            r_count += 1
            if r_count % 1e5 == 0:
                print_memory_usage("\t\t[Progress]\treads processed {} ".format(r_count))

            if read.is_unmapped or read.is_supplementary or read.is_secondary or read.is_duplicate:
                continue
            if read.is_paired and not read.is_proper_pair:
                continue
            if 'NH' in dict(read.tags) and dict(read.tags)['NH'] > 1:
                continue
    
            
            # Find the first exon that overlaps the read
            while overlapping_exon.end < read.pos and ov_index < exons.shape[0] - 2:
                ov_index += 1
                overlapping_exon = get_exon(ov_index)

            if overlapping_exon.start > read.positions[-1]:
                continue    

            ### Search for overlaps
            match_lengths, block_lengths, junction_lengths = process_CIGAR(read.cigar)
            middle_exons = range(1, len(block_lengths) - 1) 
            
            read_lengths[len(read.query_sequence)] += 1
            pos0 = read.pos

            ### Search overlap per block
            for e_i in range(len(match_lengths)): 
                match_block_s = pos0
                match_block_e = read.pos + block_lengths[e_i]
                block_gi = gi(chrom, match_block_s, match_block_e)
                
                ### Search for junction reads
                junction_gi = None
                if len(match_lengths) >= 2 and e_i <= len(match_lengths) - 2:
                    junction_s  = read.pos + block_lengths[e_i]          # 0-based, first pos on intron.  
                    junction_e  = junction_s + junction_lengths[e_i]     # 1-based, last pos on following intron. 
                    if min_intron_length <= junction_lengths[e_i] and junction_lengths[e_i] <= max_intron_length: 
                        if (match_lengths[e_i] >= overhang or e_i in middle_exons) and (match_lengths[e_i + 1] >= overhang or e_i + 1 in middle_exons):
                            junction_gi = gi(chrom, junction_s - 2, junction_e + 2)

                    pos0 = match_block_e + junction_lengths[e_i]
                
                ### Search for inclusion reads
                overlapping_exon_next = overlapping_exon
                ov_index_next = 0
                while overlapping_exon_next.start < block_gi.end:                    
                    if (match_lengths[e_i] >= overhang or e_i in middle_exons): 
                        if overlapping_exon_next.length <= micro_exon_length and block_gi.contains(overlapping_exon_next):
                            exons.at[ov_index + ov_index_next, 'INCLUSION'] += 1
                        elif overlapping_exon_next.length > micro_exon_length and block_gi.overlaps(overlapping_exon_next):
                            exons.at[ov_index + ov_index_next, 'INCLUSION'] += 1
                    ov_index_next += 1
                    overlapping_exon_next = get_exon(ov_index + ov_index_next)
                    
                    
                ### Search for exclusion reads
                overlapping_exon_next = overlapping_exon
                ov_index_next = 0
                if junction_gi is not None:
                    while overlapping_exon_next.start < junction_gi.end: 
                        if junction_gi.contains(overlapping_exon_next):               
                            exons.at[ov_index + ov_index_next, 'EXCLUSION'] += 1
                        ov_index_next += 1
                        overlapping_exon_next = get_exon(ov_index + ov_index_next)
                        
                      
        print_memory_usage("\t\t[Progress] {} {} processed.\tNumber of reads: {}".format(chrom, fstrand, r_count))
    bam.close()

    if len(read_lengths) == 0:
        raise ValueError("No reads found for this file {}".format(bam_file))

    most_used_RL = max(read_lengths, key = read_lengths.get)
    most_used_RL_prev = read_lengths[most_used_RL]/float(sum(read_lengths.values()))

    if most_used_RL_prev < 0.9:
        print_memory_usage("\tWARNING: Read length varies significantly between reads. PSI calculation may not be accurate")

    print_memory_usage("\tUsing read length = {} for PSI calculation. ({:.2%} of reads have this read length)".format(most_used_RL, most_used_RL_prev))

    return  most_used_RL



def get_PSI(exons, read_length, output_file, overhang):

    outf = open(output_file, 'w')

    outf.write("Segment_coordinates\tLength\tInclusion\tExclusion\tPSI\tPart_type\tPart_ID\n")

    for i, row in exons.iterrows():
        NIR = row["INCLUSION"]/float(read_length + row["LENGTH"] - 2*overhang)
        NER = row["EXCLUSION"]/float(read_length - 2*overhang)
        if NIR == 0 and NER == 0:
            PSI = 'nan'
        else:
            PSI = NIR/(NIR + NER)

        EXON = "{}:{}:{}:{}".format(row["CHROM"], row["START"], row["END"], row["STRAND_ORIGINAL"])

        outline = [EXON, row["LENGTH"], row["INCLUSION"], row["EXCLUSION"], PSI, row["TYPE"], row["ID"]]

        outf.write('\t'.join(map(str, outline)) + '\n')

    outf.close()



if __name__ == "__main__":
 
    description = """"PSI (Percent Spliced In) calculation of RNA-Seq data."""
    usage = """e.g. \n\t python PSI_calculator.py -b <file.bam> -a <annotation.gtf> -p <output_file>"""

    parser = OptionParser(usage = usage, description = description)
    parser.add_option("-b", "--input_bam", dest = "input_bam", 
        help = "input bam file, sorted and indexed")
    parser.add_option("-a", "--annotation_macro_exon", dest = "ma_annotation", 
        help = "annotation file (gtf format)", default = None)
    parser.add_option("-p", "--out_prefix", dest = "output_prefix", 
        help = "prefix for intermediate and output files")
    parser.add_option("--overhang", dest = "overhang", 
        help = "min length of read mapping to flanking exons in a junction [Default = 8]", default = 8)
    parser.add_option("--min_intron", dest = "min_intron_length", default = 30,
        help = "[Default = 30]")
    parser.add_option("--max_intron", dest = "max_intron_length", default = 5e5,
        help = "[Default = 500,000]")
    parser.add_option("--strandedness", dest = "strandedness", default = 'fr-secondstrand',
        help = "[Default = fr-firststrand]")
    parser.add_option("--micro_exon_length", dest = "micro_exon_length", default = 20,
        help = "Cutoff for defining micro-exons [Default = 20]")


    (options, args) = parser.parse_args()

    strandedness = options.strandedness
    min_intron_length = options.min_intron_length
    max_intron_length = options.max_intron_length
    micro_exon_length = options.micro_exon_length

    time_start = time()

    command_line = ' '.join(sys.argv)

    print_memory_usage("COMMAND:\t{}".format(command_line))
    
    print_memory_usage("Reading the annotation file {}....".format(options.ma_annotation))
    exons = read_gtf_file(options.ma_annotation, strandedness)

    print_memory_usage("Reading the bam file {}....".format(options.input_bam))
    read_length = get_coverage(exons, options.input_bam, options.overhang)

    print_memory_usage("Calculating PSI and writing to {}".format(options.output_prefix))
    get_PSI(exons, read_length, options.output_prefix, options.overhang)
	
    print_memory_usage("DONE. Run time: {:.2f} seconds.".format(time() - time_start))

