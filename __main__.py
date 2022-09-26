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
    import train_model as tm

    from sklearn.linear_model import LinearRegression
    from sklearn.neural_network import MLPRegressor
    from sklearn.svm import LinearSVR

    if args.adam == False and args.linear_regression == False and args.linearSVR == False and args.randomForest == False:
        print('Please specify at least one model type to train the fetal fraction prediction model with')
        quit()

    df = mi.combine_bins(args.folder_path,args.correlation_file, args.good_bins, args.read_lengths, args.min_cor_val, args.min_good_len, train_input = True)
    model_list = []
    model_names = []
    if args.adam == True:
        adam = MLPRegressor(hidden_layer_sizes=(5,5),solver='adam',max_iter=10000,activation="logistic",random_state=40,n_iter_no_change=5)
        model_list.append(adam)
        model_names.append("adam")

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

    tm.train_model(df, model_list, model_names, args.train_test_proportion, args.max_FFY, args.pca_components, args.print_results,args.prefix)

def predict_FF(args):
    import model_input as mi
    import predict_ff as pf
    df = mi.predict_input(args.folder_path, args.good_bins, args.read_lengths, test_input = True)
    pf.predict_FF(df,args.model_path)

def bin_coverage(args):
    import ffy
    import count_reads as cr
    contig_dictionary = cr.get_contigs(args.bam_file, args.contig_size)
    ffy.coverage(contig_dictionary,args.bam_file, args.bin_size, args.min_qual,args.max_len)

def summarise_coverage(args):
    import ffy
    ffy.summarise(args.samplesheet)

def train_ffy(args):
     import train_model as tm
     import pandas as pd
     from sklearn.linear_model import LinearRegression

     df=pd.read_csv(args.coverage_summary)
     LR = LinearRegression()
     model_list=[LR]
     model_names=['LR']

     Y=["Individual","FFY"]
     X=["Individual","FFY"]
     for col in df.columns:
         if "Y" == col[0]:
             Y.append(col)
         elif "X" == col[0]:
             X.append(col)

     df_Y=df[Y]
     df_X=df[X]

     tm.train_model(df_Y, model_list, model_names, args.train_test_proportion, args.max_FFY, args.pca_components, args.print_results,args.prefix+"_FFY")
     tm.train_model(df_X, model_list, model_names, args.train_test_proportion, args.max_FFY, args.pca_components, args.print_results,args.prefix+"_FFX")

def predict_ffy(args):
    import ffy

    header=["Individual"]
    if args.sample:
        out=[args.sample]
    else:
        out=[ args.coverage.split("/")[-1].split(".")[0] ]

    Y={}
    X={}
    n_y=1
    n_x=1

    for line in open(args.coverage):
        content=line.strip().split()
        if "Y" in content[0]:
            if len(Y) == 0:
               Y["Y"]=float(content[-1])
            else:
               Y["Y.{}".format(n_y)]=float(content[-1])
               n_y+=1

        elif "X" in content[0]:
            if len(X) == 0:
               X["X"]=float(content[-1])
            else:
               X["X.{}".format(n_x)]=float(content[-1])
               n_x+=1

    
    if args.ffy_model:
       header.append("FFY")
       pred=ffy.retrieve_ffy_model(Y,args.ffy_model)
       out.append( str(pred[0]) )

    if args.ffx_model:
       header.append("FFX")
       pred=ffy.retrieve_ffy_model(X,args.ffx_model)
       out.append( str(pred[0]) )

    print(",".join(header))
    print(",".join(out))

def main():
    if len(sys.argv) < 2:
        print("Missing submodule:")
        print("")
        print("frac_reads - count fragment length distribution")
        print("correlate_FFY - combine samples")
        print("train - train the model")
        print("predict - use the model to predict fetal fraction")
        quit()

    parser = argparse.ArgumentParser(description = "Predict fetal fraction from whole genome sequencing data.")
    
    subparsers = parser.add_subparsers()

    parser_count = subparsers.add_parser("frac_reads", help = "count help")
    parser_count.add_argument("bam_file", type=str, help = "path to bam file of individual of interest")
    parser_count.add_argument("samplesheet", type=str, help = "path to samplesheet, space separated Sample_ID and FFY")  
    parser_count.add_argument("output_folder", type=str, help = "path to folder where output files should be stored. This folder should exist prior to running the function")
    parser_count.add_argument("--bin_size","-s", type=int, default = 10000, help = "size of the chromosome bins")
    parser_count.add_argument("--min_qual","-q",  type=int, default = 40, help = "minimum mapping quality for a read to be included")
    parser_count.add_argument("--contig_size","-c",  type=int, default = 1000000, help = "minimum size of chromosomes and contigs to be included")
    parser_count.add_argument("--min_len","-i",  type=int, default = 80, help = "minimum length of reads to be included")
    parser_count.add_argument("--max_len","-a",  type=int, default = 250, help = "maximum length of reads to be included")
    parser_count.set_defaults(func=frac_reads)

    parser_correlate = subparsers.add_parser("correlate_FFY", help = "correlate help")
    parser_correlate.add_argument("folder_path", type=str, help = "string to the output folder of the count reads step")
    parser_correlate.add_argument("output_file", type=str, help = "file name that will contain the correlation data")
    parser_correlate.add_argument("--min_reads","-m",  type=int, default = 100, help = "the minimum number of reads that should be present in a bin for the bin to be considered informative")
    parser_correlate.add_argument("--min_proportion","-p",  type=float, default = 0.8, help = "the minimum proportion of times the bin should be considered informative (between different individuals) to be included")
    parser_correlate.add_argument("--min_len","-i",  type=int, default = 80, help = "minimum length of reads to be included")
    parser_correlate.add_argument("--max_len","-a",  type=int, default = 250, help = "maximum length of reads to be included")
    parser_correlate.set_defaults(func=correlate_reads)
    
    parser_input = subparsers.add_parser("train", help = "input help")
    parser_input.add_argument("correlation_file", type=str, help = "path to correlation file made in the correlate step")
    parser_input.add_argument("folder_path", type=str, help = "string to the output folder of the count reads step")
    parser_input.add_argument("good_bins", type=str, help = "txt file name that stores the good bins")
    parser_input.add_argument("read_lengths", type=str, help = "file name that stores read lengths used in model training")
    parser_input.add_argument("--min_cor_val","-c", type=float, default = 0.6, help = "minimum r value to be considered good correlation")
    parser_input.add_argument("--min_good_len","-m",  type=int, default = 5, help = "minimum number of read lengths for which a bin has good correlations for the bin to be included")
    parser_input.add_argument("--train_test_proportion","-t", type=int, default = 0.8, help = "Proportion of input used for training vs. testing")
    parser_input.add_argument("--max_FFY","-F", type=int, default = 20, help = "Max FFY of samples included in training the model. If no max is desired set to 100")
    parser_input.add_argument("--pca_components","-n", type=int, default = 10, help = "Number of PCA components data is divided into")
    parser_input.add_argument("--print_results","-p",  action="store_true", help = "Print model performance")
    parser_input.add_argument("--adam","-a",  action="store_true", help = "Use MLPRegressor in model")
    parser_input.add_argument("--linear_regression","-r",  action="store_true", help = "Use linear regressor in model")
    parser_input.add_argument("--linearSVR","-s",  action="store_true", help = "Use Linear SVR in model")
    parser_input.add_argument("--randomForest","-f",  action="store_true", help = "Use random forest in model")
    parser_input.add_argument("--prefix","-o", type=str, default = "FF", help = "Prefix of the output files")
    parser_input.set_defaults(func=train_input)

    parser_coverage = subparsers.add_parser("coverage", help = "count help")
    parser_coverage.add_argument("bam_file", type=str, help = "path to bam file of individual of interest")
    parser_coverage.add_argument("--bin_size","-s", type=int, default = 100000, help = "size of the chromosome bins")
    parser_coverage.add_argument("--min_qual","-q",  type=int, default = 40, help = "minimum mapping quality for a read to be included")
    parser_coverage.add_argument("--contig_size","-c",  type=int, default = 1000000, help = "minimum size of chromosomes and contigs to be included")
    parser_coverage.add_argument("--max_len","-a",  type=int, default = 400, help = "maximum length of reads to be included")
    parser_coverage.set_defaults(func=bin_coverage)

    parser_summarise = subparsers.add_parser("summarise_coverage", help = "correlate help")
    parser_summarise.add_argument("samplesheet", type=str, help = "path to samplesheet, space separated path and FF")  
    parser_summarise.set_defaults(func=summarise_coverage)

    parser_model_ffy = subparsers.add_parser("model_ffy", help = "model_ffy help")
    parser_model_ffy.add_argument("coverage_summary", type=str, help = "Coverage summary file produced by summarise_coverage step")  
    parser_model_ffy.add_argument("--max_FFY","-F", type=int, default = 20, help = "Max FFY of samples included in training the model. If no max is desired set to 100")
    parser_model_ffy.add_argument("--pca_components","-n", type=int, default = 4, help = "Number of PCA components data is divided into")
    parser_model_ffy.add_argument("--print_results","-p",  action="store_true", help = "Print model performance")
    parser_model_ffy.add_argument("--train_test_proportion","-t", type=float, default = 0.9, help = "Proportion of input used for training")
    parser_model_ffy.add_argument("--prefix","-o", type=str, default = "FF", help = "Prefix of the output files")
    parser_model_ffy.set_defaults(func=train_ffy)

    parser_predict_ffy = subparsers.add_parser("predict_ffy", help = "model help")
    parser_predict_ffy.add_argument("ffy_model", type=str, help = "path to FFY model file made in model_ffy step")
    parser_predict_ffy.add_argument("ffx_model", type=str, help = "path to FFY model file made in model_ffx step")
    parser_predict_ffy.add_argument("coverage", type=str, help = "txt file produced by coverage module")
    parser_predict_ffy.add_argument("--sample", type=str, help = "Sample ID")
    parser_predict_ffy.set_defaults(func=predict_ffy)

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
