#!/bin/bash
import sys
import getopt
import os
import time

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "i:")
    except getopt.GetoptError:
        print "Wrong usage, check arguments"
    for opt, arg in opts:
        if opt == "-i":
            exon_parts = arg
    start_time = time.time()
	find_introns(exon_parts)
    print "Program done"
    print "That took %i seconds" % (time.time()-start_time)


class assign(object):
    def __init__(self, lista):
        self.start, self.end = map(int, lista[3:5])
        self.part   = lista[8].split(':')[0]
        self.expart = int(lista[8].split(':')[1])


def find_introns(exon_parts):
	exon = open(exon_parts).readlines()	
	exon = [line.split() for line in exon]
	outputfile = exon_parts.replace('exonic','intronic')
	out = open(outputfile,'w')
	intronic_parts = {}
	for i in range(len(exon)-1):
		curr_part = assign(exon[i])
		next_part = assign(exon[i+1])
		chrm = exon[i][0]
		asa, strand, bsa = exon[i][5:8]
		if curr_part.part == next_part.part and next_part.start - curr_part.end > 1:
			assert next_part.expart - curr_part.expart == 1, 'Not consecutive parts'
			try: intronic_parts[curr_part.part] += 1
			except: intronic_parts[curr_part.part] = 1
			start_pos = str(curr_part.end+1)
			end_pos   = str(next_part.start-1)
			partID    = 'intron.%s:%03d' %(curr_part.part, intronic_parts[curr_part.part])
			file_name = str(sys.argv[0])
			parte   = 'intronic_part' 
			outline = [chrm, file_name, parte, start_pos, end_pos, asa, strand, bsa, partID]
			out.write('\t'.join(outline)+'\n')
	exonic.close() ; out.close()
	print 'Intronic parts made, file: %s' %(outputfile)

if __name__=='__main__':
	main(sys.argv[1:])
