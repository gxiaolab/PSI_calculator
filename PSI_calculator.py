#!/usr/bin/python
import sys
import os
import pysam 
from glob import glob 
from collections import *
from optparse import OptionParser
from os import getpid, makedirs, path, system # Kind of essential
from time import strftime, time # Not essential
from HTSeq import GenomicInterval as gi
import gc 
import itertools

import numpy as np
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
    message = strftime("[%a, %d %b %Y %H:%M:%S]\t{}\n".format(ExtraMessage))
    sys.stderr.write(message) 


def strands(read, strandedness):
    if strandedness == 'fr-firststrand':
        strand1, strand2 = '-', '+'
    elif strandedness == 'fr-secondstrand':
        strand1, strand2 = '+', '-'
    elif strandedness == 'fr-unstrand':
        strand1, strand2 = '.', '.'

    strand = ''
    ## For paired-end (PE)
    if read.is_read1:
        if not read.is_reverse:
            strand = strand1
        elif read.is_reverse:
            strand = strand2
    elif read.is_read2:
        if not read.is_reverse:
            strand = strand2
        elif read.is_reverse: 
            strand = strand1
    ## For single-end (SE)
    else:
        if not read.is_reverse:
            strand = strand1
        elif read.is_reverse:
            strand = strand2
    return strand


def process_CIGAR(cigar, pos, overhang):
    block_start = pos
    all_blocks = [] 
    read_match   = 0 
    genome_match = pos
    for (ctype, length) in cigar:
        if ctype == 0: #M
            read_match   += length
            genome_match += length
        elif ctype == 3: #N
            all_blocks.append((block_start, genome_match, read_match))
            genome_match += length
            block_start = genome_match
            read_match = 0
        elif ctype == 2: #D
            genome_match += length
    all_blocks.append((block_start, genome_match, read_match))

    if all_blocks and all_blocks[0][2] < overhang:
        all_blocks = all_blocks[1:]

    if all_blocks and all_blocks[-1][2] < overhang:
        all_blocks = all_blocks[:-1]
 
    return all_blocks


START_INDEX = 0
END_INDEX = 1
ID_INDEX = 2
TYPE_INDEX = 3
LENGTH_INDEX = 4
STRAND_ORIGINAL_INDEX = 5
INCLUSION_INDEX = 6
EXCLUSION_INDEX = 7


def read_gtf_file(gtf_file, strandedness):
    exons = defaultdict(list)
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

            exons[(chrom, strand)].append([feature_start - 1, feature_end, feature_id, feature_type, feature_length, original_strand, 0, 0])

    for contig, coord_list in exons.items():
        exons[contig] = sorted(coord_list, key = lambda x: (x[0], x[1]))

    return exons

      


def get_coverage(exons, bam_file,  overhang):

    read_lengths = defaultdict(int)

    references_not_found = set()

    for contig in exons:

        (chrom, fstrand) = contig

        last = len(exons[contig])

        def get_exon(overlap_index):
            if overlap_index >= last:
                return gi(chrom, sys.maxint - 2, sys.maxint)
            else:
                me_start = exons[contig][overlap_index][START_INDEX]
                me_end   = exons[contig][overlap_index][END_INDEX]
                return gi(chrom, me_start, me_end)

        r_count = 0
        first_read_index = 0
        overlapping_exon = get_exon(first_read_index)

        bam = pysam.Samfile(bam_file, 'rb')
        
        if chrom not in bam.references:
            references_not_found.add(chrom)
            bam.close()
            continue  

        for read in bam.fetch(chrom, until_eof = True): 
            
            strand = strands(read, strandedness)

            if strand != fstrand:
                continue

            r_count += 1
            if r_count % int(5e5) == 0:
                print_memory_usage("[Progress]\treference {0} reads processed {1:,}...".format(chrom, r_count))
                

            if read.is_unmapped or read.is_supplementary or read.is_secondary or read.is_duplicate:
                continue
            if read.is_paired and not read.is_proper_pair:
                continue
            if 'NH' in dict(read.tags) and dict(read.tags)['NH'] > 1:
                continue
    
            read_lengths[len(read.query_sequence)] += 1

            # Find the first exon that overlaps the read
            while overlapping_exon.end < read.pos:
                first_read_index += 1
                overlapping_exon = get_exon(first_read_index)

            if overlapping_exon.start > read.positions[-1]:
                continue    

            ### Search for overlaps
            all_blocks = process_CIGAR(read.cigar, read.pos, overhang)
            
            middle_exons = range(1, len(all_blocks) - 1)   

            ### Search overlap per block
            for e_i in range(len(all_blocks)): 

                block_s, block_e, block_match_length = all_blocks[e_i]

                block_gi = gi(chrom, block_s, block_e)

                ### Search for junction reads
                junction_gi = None

                if len(all_blocks) >= 2 and e_i <= len(all_blocks) - 2:
                    junction_s  = all_blocks[e_i][1]       # 0-based, first pos on intron.  
                    junction_e  = all_blocks[e_i + 1][0]   # 1-based, last pos on following intron. 
                    junction_length = junction_e - junction_s

                    if min_intron_length <= junction_length and junction_length <= max_intron_length:     
                        junction_gi = gi(chrom, junction_s - 1, junction_e + 1)
                
                ### Search for exons that overlap the read block
                first_block_index = first_read_index

                overlapping_exon_next = get_exon(first_block_index) #overlapping_exon
           
                while overlapping_exon_next.start < block_gi.end:  
                       
                    if overlapping_exon_next.length <= micro_exon_length and block_gi.contains(overlapping_exon_next) and e_i in middle_exons:
                        exons[contig][first_block_index][INCLUSION_INDEX] += 1
                    if overlapping_exon_next.length > micro_exon_length and block_gi.overlaps(overlapping_exon_next):
                        overlap_length = min(block_gi.end - overlapping_exon_next.start, overlapping_exon_next.end - block_gi.start)
                        if overlap_length >= overhang:
                            exons[contig][first_block_index][INCLUSION_INDEX] += 1
                    
                    first_block_index += 1 
                    overlapping_exon_next = get_exon(first_block_index)

                ### Search for exclusion reads
                if junction_gi is not None:
                    while overlapping_exon_next.start < junction_gi.end: 
                        if junction_gi.contains(overlapping_exon_next):               
                            exons[contig][first_block_index][EXCLUSION_INDEX] += 1
                        first_block_index += 1 
                        overlapping_exon_next = get_exon(first_block_index)

        bam.close()     
        gc.collect()       
        print_memory_usage("[Progress] {0} {1} done.\tTotal number of reads: {2:,}.".format(chrom, fstrand, r_count))
   
    print_memory_usage("\t\tWARNING: references not found in bam file: {}".format(', '.join(references_not_found)))


    if len(read_lengths) == 0:
        raise ValueError("No reads found for this file {}".format(bam_file))

    most_used_RL = max(read_lengths, key = read_lengths.get)
    most_used_RL_prev = read_lengths[most_used_RL]/float(sum(read_lengths.values()))

    if most_used_RL_prev < 0.9:
        print_memory_usage("\tWARNING: Read length varies significantly between reads. PSI calculation may not be accurate")

    print_memory_usage("\tUsing read length = {} for PSI calculation. ({:.2%} of reads have this read length)".format(most_used_RL, most_used_RL_prev))

    return most_used_RL



def get_PSI(exons, read_length, output_file, overhang):

    outf = open(output_file, 'w')

    outf.write("Segment_coordinates\tLength\tInclusion\tExclusion\tPSI\tPart_type\tPart_ID\n")

    all_outlines = []
    for (chrom, strand), coord_list in exons.items():
        for row in coord_list:
            if row[LENGTH_INDEX] > micro_exon_length:
                NIR = row[INCLUSION_INDEX]/float(read_length + row[LENGTH_INDEX] - 2*overhang)
                NER = row[EXCLUSION_INDEX]/float(read_length - 2*overhang)
            else:
                NIR = row[INCLUSION_INDEX]/float(read_length + row[LENGTH_INDEX] - 1)
                NER = row[EXCLUSION_INDEX]/float(read_length - 1)

            if NIR == 0 and NER == 0:
                PSI = 'nan'
            else:
                PSI = NIR/(NIR + NER)

            EXON = "{}:{}:{}:{}".format(chrom, row[START_INDEX], row[END_INDEX], row[STRAND_ORIGINAL_INDEX])

            outline = [chrom, row[START_INDEX], EXON, row[LENGTH_INDEX], row[INCLUSION_INDEX], row[EXCLUSION_INDEX], PSI, row[TYPE_INDEX], row[ID_INDEX]]
            all_outlines.append(outline)

    all_outlines.sort(key = lambda x: (x[0], x[1]))
    for outline in all_outlines:
        outf.write('\t'.join(map(str, outline[2:])) + '\n')
    outf.close()



if __name__ == "__main__":
 
    description = """"PSI (Percent Spliced In) calculation of RNA-Seq data."""
    usage = """e.g. \n\t python PSI_calculator.py -b <file.bam> -a <annotation.gtf> -p <output_file>"""

    parser = OptionParser(usage = usage, description = description)
    parser.add_option("-b", "--input_bam", dest = "input_bam", 
        help = "input bam file, sorted and indexed")
    parser.add_option("-a", "--annotation", dest = "annotation", 
        help = "annotation file (gtf format)", default = None)
    parser.add_option("-p", "--out_prefix", dest = "output_prefix", 
        help = "prefix for intermediate and output files")
    parser.add_option("--overhang", dest = "overhang", 
        help = "min length of read mapping to flanking exons in a junction [Default = 8]", default = 8)
    parser.add_option("--min_intron", dest = "min_intron_length", default = 30,
        help = "[Default = 30]")
    parser.add_option("--max_intron", dest = "max_intron_length", default = 5e5,
        help = "[Default = 500,000]")
    parser.add_option("--strandedness", dest = "strandedness", default = 'fr-firststrand',
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
    
    print_memory_usage("Reading the annotation file {}....".format(options.annotation))
    exons = read_gtf_file(options.annotation, strandedness)

    print_memory_usage("Reading the bam file {}....".format(options.input_bam))
    read_length = get_coverage(exons, options.input_bam, options.overhang)

    print_memory_usage("Calculating PSI and writing to {}".format(options.output_prefix))
    get_PSI(exons, read_length, options.output_prefix, options.overhang)
	
    print_memory_usage("DONE. Run time: {:.2f} seconds.".format(time() - time_start))

