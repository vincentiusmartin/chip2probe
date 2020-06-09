# chip2probe

## chip2probe available modules
1. probe_generator
2. training_gen
3. modeler
4. bio

## imads directory is for generating imads model

## README

TODO
Use Kompas to generate site_all.tsv without imads score

1. Download ChiP-seq files
Background

We generate probes using sequences from in vivo data (i.e. ChIP-seq) to enable binding measurement of the sequences that TFs do encounter in the cells. Therefore, the first step in the pipeline is to choose which ChIP-seq data to use.

Data sources

There are some ChIP-seq data sources available:
1. Encode: https://www.encodeproject.org/
2. Cistrome: http://cistrome.org/db/#/
3. ReMap: http://tagc.univ-mrs.fr/remap/

Pipeline requirement

For the data download part, the program is able to facilitate downloads by enabling user to simply input the urls. The subsequent process will be tailored based on each of the data sources inputted by the user. Look at dictinput_example.py for input examples for the data downloading process. Create a similar file with user inputs in a file called 'dictinput.py'

2. Quality Check

At this point, the data is usually already in a narrow peak format either as the result of peak calling or the raw data itself is already in this format. The next step is to check the quality of the data, for example: if the ChIP-seq data has a good quality, then it should have a good correlation between replicates. Also we want to do some information such as: how many predicted binding sites can we find in each peak? All of this information will help us in determining which data to use also later to go back to the data after we do some analyses to our data.

Quality Check Output
By default the program will provide all information that it could gather from the data.

Requirements
I. Input parser:

(these two functions assume the columns format, make it general?)

1. read.narrow.peak:

2. read.pileup:

Running

In run_all.py,
In run_all.sh, edit outdir to be the path to the directory the ChIP-seq files are to be downloaded to. Index will be the index of the dictionary in inList list defined in dictinput.py.

Output
sites_all.tsv

3. Probe filter
analysis_path -> path to site_all.tsv

Download libsvm:
https://www.csie.ntu.edu.tw/~cjlin/libsvm/oldfiles/index-1.0.html
libsvm >= 3.24
Unzip and rename to 'libsvm'
Put libsvm under sitespredict directory
Set up libsvm:
Run 'make' in libsvm directory
Run 'make' in libsvm/python directory
