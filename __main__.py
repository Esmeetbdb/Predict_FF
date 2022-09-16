import argparse
import sys
import os

def frac_reads(args):
    import count_reads as cr
    contig_dictionary = cr.get_contigs(args.bam_file, args.contig_size)
    counted_reads = cr.extract_read_length(args.bam_file,
                        contig_dictionary,
                        args.bin_size,
                        args.min_qual,
                        args.min_len,
                        args.max_len)
    ffy_dict = cr.fetal_fraction_dict(args.samplesheet)
    cr.write_output(args.bam_file,
                 ffy_dict,
                 counted_reads,
                 args.output_folder,
                 args.min_len,
                 args.max_len)

def correlate_reads(args):
    import correlate_FFY as cF
    corr_dict = cF.get_correlations(args.folder_path, args.min_len, args.max_len, args.min_reads, args.min_proportion)
    cF.correlation_csv(corr_dict, args.min_len, args.max_len, args.output_file)

def train_input(args):
    import model_input as mi
    df = combine_bins(args.folder_path,args.correlation_file, args.good_bins, args.read_lengths, args.min_cor_val, args.min_good_len, train_input = True)
    model_list = []
    model_names = []
    if args.adam == True:
        adam = MLPRegressor(hidden_layer_sizes=(5,5),solver='sgd',max_iter=10000,activation="logistic",random_state=40,n_iter_no_change=5)
        model_list.append(adam)
        model_names.appen("adam")

    if args.linear_regression == True:
        LR = LinearRegression()
        model_list.append(LR)
        model_names.append('LR')

    if args.linearSVR == True:
        lvsm = LinearSVR(random_state=0, tol=1e-5,C=10,max_iter=1000000)
        model_list.append(lvsm)
        model_names.append('lvsm')

    if args.randomForest == True:
        randomforest=RandomForestRegressor(max_depth=10, random_state=0,n_estimators=100)
        model_list.append(randomforest)
        model_names.append('randomforest')

    elif args.adam == False and args.linear_regression == False and args.linearSVR == False and args.randomForest == False:
        print('Please specify at least one model type to train the fetal fraction prediction model with')

    train_model(df, model_list, model_names, args.train_test_proportion, args.max_FFY, args.pca_components, args.print_results)

def predict_FF(args):
    import model_input as mi
    import predict_ff as pf
    df = mi.predict_input(args.folder_path, args.good_bins, args.read_lengths, test_input = True)
    pf.predict_FF(df)

def main():
    parser = argparse.ArgumentParser(description = "Predict fetal fraction from whole genome sequencing data.")
    
    subparsers = parser.add_subparsers()

    parser_count = subparsers.add_parser("frac_reads", help = "count help")
    parser_count.add_argument("bam_file", type=str, help = "path to bam file of individual of interest")
    parser_count.add_argument("samplesheet", type=str, help = "path to samplesheet, space separated Sample_ID and FFY")  
    parser_count.add_argument("output_folder", type=str, help = "path to folder where output files should be stored. This folder should exist prior to running the function")
    parser_count.add_argument("--bin_size","-s", type=int, default = 10000, help = "size of the chromosome bins")
    parser_count.add_argument("--min_qual","-q",  type=int, default = 10, help = "minimum mapping quality for a read to be included")
    parser_count.add_argument("--contig_size","-c",  type=int, default = 1000000, help = "minimum size of chromosomes and contigs to be included")
    parser_count.add_argument("--min_len","-i",  type=int, default = 80, help = "minimum length of reads to be included")
    parser_count.add_argument("--max_len","-a",  type=int, default = 250, help = "maximum length of reads to be included")
    parser_count.set_defaults(func=frac_reads)

    parser_correlate = subparsers.add_parser("correlate_FFY", help = "correlate help")
    parser_correlate.add_argument("folder_path", type=str, help = "string to the output folder of the count reads step")
    parser_correlate.add_argument("output_file", type=str, help = "file name that will contain the correlation data")
    parser_count.add_argument("--min_reads","-m",  type=int, default = 100, help = "the minimum number of reads that should be present in a bin for the bin to be considered informative")
    parser_count.add_argument("--min_proportion","-p",  type=float, default = 0.8, help = "the minimum proportion of times the bin should be considered informative (between different individuals) to be included")
    parser_correlate.add_argument("--min_len","-i",  type=int, default = 80, help = "minimum length of reads to be included")
    parser_correlate.add_argument("--max_len","-a",  type=int, default = 250, help = "maximum length of reads to be included")
    parser_correlate.set_defaults(func=correlate_reads)
    
    parser_input = subparsers.add_parser("train", help = "input help")
    parser_input.add_argument("correlation_file", type=str, help = "path to correlation file made in the correlate step")
    parser_input.add_argument("folder_path", type=str, help = "string to the output folder of the count reads step")
    parser_input.add_argument("good_bins", type=str, help = "txt file name that stores the good bins")
    parser_input.add_argument("read_lengths", type=str, help = "file name that stores read lengths used in model training")
    parser_input.add_argument("--min_cor_val","-c", type=int, default = 0.6, help = "minimum r value to be considered good correlation")
    parser_input.add_argument("--min_good_len","-m",  type=int, default = 5, help = "minimum number of read lengths for which a bin has good correlations for the bin to be included")
    parser_input.add_argument("--train_test_proportion","-t", type=int, default = 0.8, help = "Proportion of input used for training vs. testing")
    parser_input.add_argument("--max_FFY","-F", type=int, default = 20, help = "Max FFY of samples included in training the model. If no max is desired set to 100")
    parser_input.add_argument("--pca_components","-n", type=int, default = 10, help = "Number of PCA components data is divided into")
    parser_input.add_argument("--print","-p",  action="store_true", help = "Print model performance")
    parser_input.add_argument("--adam","-a",  action="store_true", help = "Use MLPRegressor in model")
    parser_input.add_argument("--linear_regression","-r",  action="store_true", help = "Use linear regressor in model")
    parser_input.add_argument("--linearSVR","-s",  action="store_true", help = "Use Linear SVR in model")
    parser_input.add_argument("--randomForest","-o",  action="store_true", help = "Use random forest in model")
    parser_input.set_defaults(func=train_input)

    parser_predict = subparsers.add_parser("predict", help = "model help")
    parser_predict.add_argument("model_path", type=str, help = "path to model file made in train_model step")
    parser_predict.add_argument("folder_path", type=str, help = "string to the output folder of the count reads step")
    parser_predict.add_argument("good_bins", type=str, help = "txt file name that stores the good bins")
    parser_predict.add_argument("read_lengths", type=str, help = "file name that stores read lengths used in model training")
    parser_predict.set_defaults(func=predict_FF)

    args = parser.parse_args(sys.argv[1:])
    args.func(args)
if __name__ == '__main__':
    main()
