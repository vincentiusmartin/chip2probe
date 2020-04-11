# chip2probe

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

Running

In run_all.py, 
In run_all.sh, edit outdir to be the path to the directory the ChIP-seq files are to be downloaded to. Index will be the index of the dictionary in inList list defined in dictinput.py.

Output
sites_all.tsv

2. 
