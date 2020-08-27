# iMADS model generator

## Input paramater for inputdict.py
1. pbmdata: path to the pbm file
2. column_train: which column has the intensity value (right now we assume the input has been normalized)
3. kmers: which k to use
4. width: imdas model width
5. corelist: list of the cores for the corresponding tf
6. tfname: name of the transcription factor, e.g. "Runx1"
7. grid: SVR parameter grid
8. numfold: fold for cross validation
9. corepos: position of the core (center, left, right)
10. outdir: output folder
11. logit: do logistic transformation on the input (True/False)

## Running locally
1. Edit inputdict.py
2. `python imads_main.py`

## Running in HARDAC
1. Edit inputdict.py
2. On SLURM based environment, run: `./run_all.sh`

## Output
The script will output an SVR model file which can be loaded by `svmutil.svm_load_model`.
The model can also be used for sitespredict, see `chip2probe/main/sitespredict/plotseq.py`.
