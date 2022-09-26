import pysam
import os
import pandas as pd
import numpy as np

def get_contigs(bam_file, min_contig_size=1000000):
    ''' Get all chromosomes that should be included in model building and fetal fraction prediction'''

    sam_file = pysam.AlignmentFile(bam_file, "rb")
    bam_header = sam_file.header
    contig_length = {}

    for contig in bam_header["SQ"]:
        if int(contig["LN"]) > min_contig_size:
            contig_length[contig["SN"]] = contig["LN"]
    return contig_length

def read_ok(read, min_mapping_quality):
    ''' Check if the read has a high enough quality and if botht he read and it's mate are mapped.'''

    if read.is_read1 and not read.mate_is_unmapped and not read.is_unmapped and not read.is_duplicate:
        if read.mapping_quality > min_mapping_quality and ( (read.is_reverse and not read.mate_is_reverse) or (not read.is_reverse and read.mate_is_reverse) ):
            return True

def extract_read_length(bam_file,
                        contig_dictionary,
                        chromosome_bin_size,
                        min_mapping_quality,
                        min_read_length,
                        max_read_length):
    ''' Counts the prevalence of all read lengths per chromosome bin of specified size. Smaller bins will take longer but may result in more precise results. Too small bins
        results in no correlations. '''

    samfile = pysam.AlignmentFile(bam_file, "r")
    insert_size_per_region = {}
    frac_per_region = {}

    for contig in contig_dictionary:
        for i in range(0,contig_dictionary[contig],chromosome_bin_size):
            idx=contig+'_'+str(i)+'_'+str(i+chromosome_bin_size)
            insert_size_per_region[idx] = {}
            frac_per_region[idx] = {}
            for x in range(min_read_length,max_read_length+1):
                insert_size_per_region[idx][x] = 0
                frac_per_region[idx][x] = 0.0
            total_len = 0
            for read in samfile.fetch(contig, i, (i+chromosome_bin_size)):
                read_length = read.template_length
                if read_length >= min_read_length and read_ok(read, min_mapping_quality) and read_length <= max_read_length:
                    total_len += 1
                    insert_size_per_region[idx][read_length] += 1
            frac_per_region[idx]['total_reads'] = total_len
            if total_len == 0:
                continue
            for j in insert_size_per_region[idx]:
                frac_per_region[idx][j] = insert_size_per_region[idx][j] / float(total_len)
    return frac_per_region

def fetal_fraction_dict(samplesheet):
    ''' Get FFY from the samplesheet '''

    fetal_fraction = {}
    for line in open(samplesheet):
        content = line.strip().split(' ')
        ind = content[0]
        ffy = content[1]
        fetal_fraction[ind] = ffy
    return fetal_fraction

def write_output(bam_file,
                 FFY_dictionary,
                 frac_all_dictionary,
                 output_folder,
                 min_read_length,
                 max_read_length):
    ''' Write the output file that contains the read length counts per chromosome bin, the total reads per bin, FFY and sample ID '''

    temp = bam_file.split('/')
    ind = temp[-1].replace('.bam', '')

    file = open(output_folder+'/'+ind+'.bamoutput.csv', 'w')
    file.writelines('Individual,FFY,bin,total_reads')

    for i in range(min_read_length, max_read_length):
        file.writelines(','+str(i))
    file.writelines('\n')
    
    if ind in FFY_dictionary:
        ffy = FFY_dictionary[ind]
        for bins in frac_all_dictionary:
            file.writelines(ind+','+ffy + ',' + bins + ',' + str(frac_all_dictionary[bins]['total_reads']))
            for i in range(min_read_length,max_read_length):
                if i not in frac_all_dictionary[bins]:
                    file.writelines(',' + '0.0')
                else:
                    file.writelines(',' + str(frac_all_dictionary[bins][i]))
            file.writelines('\n')
    file.close()
