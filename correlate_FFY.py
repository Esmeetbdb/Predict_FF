import os
import pandas as pd
import numpy as np

def file_is_male(file_name):
    ''' Esure all included samples are male '''
    for line in open(file_name):
        if 'vidual' in line:
            continue
        content = line.strip().split(',')
        ffy = content[1]
        if ffy == '0.0':
            return None
        else:
            return ffy

def exclude_bins(file_names_list, min_reads, min_proportion):
    ''' exlude bins that do not have enough reads per bin in a large enough fraction of bins '''
    count_bins = {}
    good_bins = []
    for file in file_names:
        for line in open(file):
            if 'vidual' in line:
                continue
            content = line.strip().split(',')
            bins = content[2]
            total_reads = content[3]
            if bins not in count_bins:
                count_bins[bins] = 0
            if total_reads >= min_reads:
                count_bins[bins] +=1 
    for bins in count_bins:
        if count_bins[bins]/float(len(file_names_list)) >= min_proportion:
            good_bins.append(bins)
    return good_bins

def get_files(path_to_csv, train = False, test = False):
    ''' Get only files belonging to male samples '''
    paths = os.listdir(path_to_csv)
    file_names = []
    if train == True:
        for file in paths:
            files = path_to_csv+'/'+file
            if files.endswith('.bamoutput.csv') and os.path.isfile(files):
                male_ffy = file_is_male(files)
                if male_ffy != None:
                    file_names.append(files)

    elif test ==True:
        for file in paths:
            files = path_to_csv+'/'+file
            if files.endswith('.bamoutput.csv') and os.path.isfile(files):
                file_names.append(files)
    return file_names

def get_ffy(file_names_list):
    ''' Save ffy of each sample in an ordered list '''
    ffy_list = []
    for file in file_names_list:
        ffy = file_is_male(file)
        ffy_list.append(ffy)
    return ffy_list

def get_frac_per_bin(bin_list, file_names_list):
    frac_per_len = {}
    for chrom_bin in bin_list:
        frac_per_len[chrom_bin] = {}
    for line in open(file_names_list[0]):
        if 'vidual' in line:
            content = line.strip().split(',')
            read_lengths = content[4:]
            continue

    for file in file_names_list:
        for line in open(file):
            if 'vidual' in line:
                continue
            content = line.strip().split(',')
            chrom_bin = content[2]
            if chrom_bin not in bin_list:
                continue
            for i in range(4,4+len(read_lengths)):
                if read_lengths[i-4] not in frac_per_len[chrom_bin]:
                    frac_per_len[chrom_bin][read_lengths[i-4]] = [content[i]]
                else:
                    frac_per_len[chrom_bin][read_lengths[i-4]].append(content[i])

    return frac_per_len

def get_correlations(path_to_csv, min_read_len, max_read_len, min_reads, min_proportion):
    ''' Get correlations per read length per bin of frac read length to FFY '''
    corr_per_len = {}

    file_names = get_files(path_to_csv, train = True)
    ffy_list = get_ffy(file_names)
    bin_list = exclude_bins(file_names,min_reads, min_proportion)
    frac_per_bin = get_frac_per_bin(bin_list, file_names)

    var1 = np.zeros(len(file_names))
    for i in range(len(var1)):
        var1[i]=ffy_list[i]

    var2 = np.zeros(len(file_names))
    for chrom_bin in frac_per_len:
        corr_per_len[chrom_bin] = {}
        for i in range(min_read_len,max_read_len):
            frac_list = frac_per_len[chrom_bin][i]

            for j in range(len(var2)):
                var2[j] = frac_list[j]
            R_number = np.corrcoef(var1, var2)[0,1]
            corr_per_len[chrom_bin][i] = R_number


    return corr_per_len

def correlation_csv(corr_dict, min_read_length, max_read_length, output_name):
    file = open(output_name, 'w')
    for i in range(min_read_length, max_read_length):
        file.writelines(',' + str(i))
    file.writelines('\n')

    for chrom_bin in corr_dict:
        file.writelines(chrom_bin)
        for i in range(min_read_length, max_read_length):
            if i in corr_dict[chrom_bin]:
                file.writelines(',' + str(corr_dict[chrom_bin][i]))
            else:
                file.writelines(', -')
        file.writelines('\n')


def find_best_correlations(correlation_dict):
    best_cor = {}
    for chrom_bin in correlation_dict:
        best_cor[chrom_bin] = {}
        for read_len in correlation_dict[chrom_bin]:
            R_number = correlation_dict[chrom_bin][read_len]
            if R_number >= 0.7 or R_number <= -0.7 :
                best_cor[chrom_bin][read_len] = R_number
    return best_cor
