# iMADS model generator

## Input paramater for inputdict.py
1. pbmdata: path to the pbm file
2. column_train: which column has the intensity value (right now we assume the input has been normalized)
3. kmers: which k to use
4. width: imads model width, usually 20 (not the the sequence length)
5. corelist: list of the cores for the corresponding tf
6. tfname: name of the transcription factor, e.g. "Runx1"
7. grid: SVR parameter grid
8. numfold: fold for cross validation
9. corepos: position of the core (center, left, right)
10. outdir: output folder
11. logit: do logistic transformation on the input (True/False)

# Preparing input
1. Rename inputdict_example.py to `inputdict.py`
2. Edit `inputdict.py`

## Running locally
`python imads_main.py`

## Running in HARDAC
1. Need to set things up, using interactive session: `srun -p interactive --mem 32GB --pty /bin/bash`
2. `module load Anaconda3/4.3.0-gcb01`
3. `module load gcc/4.7.4-fasrc01`
4. Install libsvm: `conda install libsvm`
5. Install other packages: `conda install [package_name]`
6. Running in SLURM: `./run_all.sh`

## Output
The script will output an SVR model file which can be loaded by `svmutil.svm_load_model`.
The model can also be used for sitespredict, see `chip2probe/main/sitespredict/plotseq.py`.
