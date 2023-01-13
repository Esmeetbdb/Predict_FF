# Predict_FF
Predict fetal fraction from sequenced DNA  extracted from maternal blood

## Run
### Step 1: Counting read length fractions
Before starting model training or prediction, BAM files need to be converted to CSV files containing read length fractions. This can be done using:

```
python __main__.py frac_reads bam_file samplesheet output_folder [--optional arguments]
```

Optional parameters are:
Optional argument | Function | Default
---|---|---
--bin_size x | Change the size of the chromosome bins | 10000
--min_qual | Minimum mapping quality for a read to be included in further analysis | 10
--contig_size | Minimum size of chromosome and contigs to be included in mapping | 1000000
--min_len | Minimum length of reads to be included | 80
--max_len | Maximum length of reads to be included | 250

For more information about running this command use
 ```
 python __main__.py frac_reads --help 
 ```

### Step 2: Correlate FFY to read length fractions
FFY is then correlated to read length to find regions that can predict FFY and to remove regions that do not contain useful information

```
python __main__.py correalte_FFY folder_path output_file
```

Optional parameters are:
Optional arhument | Function | Default
---|---|---
--min_reads | the minimum number of reads that should be present in a bin for the bin to be considered informative | 100
--min_proportion | the minimum proportion of times the bin should be considered informative (between different individuals) to be included | 0.8
--min_len | Minimum length of reads to be included | 80
--max_len | Maximum length of reads to be included | 250

For more information about running this command use
 ```
 python __main__.py correalte_FFY --help 
 ```
 
 ### Step 3: Training the model
 The model is then trained using FFY data and the data generated in the previous steps.
 
 ```
python __main__.py train correlation_file folder_path good_bins read_lengths
```
Optional parameters are:
Optional arhument | Function | Default
---|---|---
--min_cor_val | minimum r value to be considered good correlation | 0.6
--min_good_len | minimum number of read lengths for which a bin has good correlations for the bin to be included | 5
--train_test_proportion | Proportion of input used for training vs. testing | 0.8
--max_FFY | Max FFY of samples included in training the model. If no max is desired set to 100 | 20
--pca_components | "Number of PCA components data is divided into | 10
--print | Print model performance | True
--adam | Use MLPRegressor in model | True
--linear_regression | Use linear regressor in model | True
--linearSVR | Use Linear SVR in model | True
--randomForest | Use random forest in model | True

For more information about running this command use
 ```
 python __main__.py train --help 
```
### Step 4: Use the model to predict fetal fraction
The model trained in step 3 can now be used to predict fetal fraction of future samples
```
python __main__.py predict model_path folder_path good_bins read_lengths
```
For more information about running this command use
 ```
 python __main__.py predict --help 
```
