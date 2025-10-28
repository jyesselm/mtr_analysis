# mtr_analysis

[![PYPI status]( https://badge.fury.io/py/mtr_analysis.png)](http://badge.fury.io/py/mtr_analysis)[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

some simple analysis scripts for processing mtr ribozyme data

## Install

To install mtr_analysis 

```shell
git clone https://github.com/jyesselm/mtr_analysis
cd mtr_analysis
conda env create -f environment.yml
conda activate mtr-analysis
python -m pip install -e .

# can test commands in example directory

```
## How to use 

### run rna-map to get bitvector files for each construct
This is the longest step will take a while to run > 1 hour.
```shell 
# download the demultiplexed data and data.csv from $NRDSTOR/run_name 
mtr-analysis run-rna-map data.csv
# expected output:

                construct  barcode_seq  time
  mtr1_mut_lib_t15min_II_ CTGCGTGCAAAC    15
  mtr1_mut_lib_t60min_II_ GCAAATGTGCTA    60
 mtr1_mut_lib_t180min_II_ AAGGACCACTGG   180
 mtr1_mut_lib_t420min_II_ CGGGCACGGCGG   420
mtr1_mut_lib_t1440min_II_ CTAGCAATGTGA  1440
# in between here you will see the output of the rna-map command, and the summary.csv file will be saved in data/construct_name
# at the end you will see the summary.csv file for all the constructs
                   construct  time    reads  aligned
0    mtr1_mut_lib_t15min_II_    15  8372027    98.38
1    mtr1_mut_lib_t60min_II_    60  7073223    98.88
2   mtr1_mut_lib_t180min_II_   180  5069385    98.63
3   mtr1_mut_lib_t420min_II_   420  9167912    98.58
4  mtr1_mut_lib_t1440min_II_  1440      654    97.25

```
### get the mutation fractions for each construct
```shell
mtr-analysis get-mutation-fractions
# this currently doesnt have any output but will create the mut_fractions.csv file in each construct directory 
# it will also create the all_mut_fractions.csv file in the current directory which 

```

### fit the mutation fractions to a monoexponential model
This will fit the mutation fractions to a monoexponential model and plot the results.

It will create the mut_kinetics.csv file in the current directory which will have the y_max, k, and k_std for each mutation.

It will also create the plots/mut.png file in the current directory which will have the plot of the mutation fractions and the fit.
```shell
mtr-analysis fit-mut-fractions
```

