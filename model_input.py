import os
import pandas as pd
from operator import itemgetter
import random as rand
import numpy
from correlate_FFY import file_is_male
from correlate_FFY import get_files


def get_good_bins(correlation_vals_file, min_correlation_val, min_good_lengths):
    good_bins = []
    for line in open(correlation_vals_file):
        if '80,81' in line:
            continue
        content = line.strip().split(',')
        bins = content[0]
        n = 0
        for item in content[1:]:
            if abs(float(item)) >= min_correlation_val:
                n += 1
        if n >= min_good_lengths:
            good_bins.append(bins)
    return good_bins

def combine_bins(frac_file_dir,correlation_vals_file, good_bin_csv, read_length_csv, min_correlation_val, min_good_lengths, test_input = False, train_input = False):

    if train_input == True:
        files = get_files(frac_file_dir, train = True)

        good_bins = get_good_bins(correlation_vals_file, min_correlation_val, min_good_lengths)
        with open(good_bin_csv, 'w') as file:
            for item in good_bins:
                file.writelines(item+'\n')

        for line in open(files[0]):
            if 'vidual' in line:
                content = line.strip().split(',')
                read_lengths = content[3:]
                break
        with open(read_length_csv, 'w') as file:
            for item in read_lengths:
                file.writelines(item + '\n')
    elif test_input == True:
        files = get_files(frac_file_dir, test = True)
        good_bins = []
        read_lengths = []
        for line in open(good_bin_csv):
            good_bins.append(line.strip())
        for line in open(read_length_csv):
            read_lengths.append(line.strip())

    avg_frac = {}
       
    for file in files:
        ind = str(file.split('/')[-1]).replace('.bamoutput.csv','')
        avg_frac[ind] = {}
        df = pd.read_csv(file)
        df_good = df[df['bin'].isin(good_bins)]
        FFY = df_good['FFY'].mean()
        avg_frac[ind]['FFY'] = FFY
        for rl in read_lengths:
            avg = df_good[rl].mean()
            avg_frac[ind][rl] = avg
    df = pd.DataFrame.from_dict(avg_frac, orient = 'columns')
    return df

def predict_input(frac_file_dir, good_bin_csv, read_length_csv, test_input = True):
    combined_bins = combine_bins(frac_file_dir, None, good_bin_csv, read_length_csv, None, None, test_input = True)
    df = pd.DataFrame.from_dict(combined_bins, orient='columns')
    return df


def combined_bins_output(frac_file_dir,correlation_vals_file, good_bins_csv, read_length_csv, output_file, min_correlation_val, min_good_lengths, test_input = False, train_input = False):
    combined_bins = combine_bins(frac_file_dir,correlation_vals_file, good_bin_csv, read_length_csv,min_correlation_val, min_good_lengths, test_input = test_input, train_input = train_input)
    read_lengths = []
    file = open(output_file, 'w')
    inds = list(combined_bins.keys())
    file.writelines('Individual,FFY')
    for rl in combined_bins[inds[0]]:
        if rl == 'FFY':
            continue
        file.writelines(',' + str(rl))
        read_lengths.append(rl)
    file.writelines('\n')
    print(len(read_lengths))
    for ind in inds:
        FFY = combined_bins[ind]['FFY']
        file.writelines(ind + ',' + str(FFY))
        for rl in read_lengths:
            file.writelines(',' + str(combined_bins[ind][rl]))
        file.writelines('\n')


