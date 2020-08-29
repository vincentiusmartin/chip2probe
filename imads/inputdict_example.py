param = {
        "pbmdata" : "/path/to/input.txt",
        "column_train" : "Ets1_100nM",
        "kmers" : [1,2,3],
        "width" : 12,
        "corelist" : ["GGAA", "GGAT"],
        "tfname" : "Ets1",
        "grid" : { #SVR feature
            "c": [0.01, 0.05, 0.1, 0.5, 1, 5, 10], # cost
            "g": [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 1, 10], # gamma for epsilon SVR
            "p": [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1], # epsilon, linear don't use: 0.01
            "t": [2] # kernel type: 0:linear, 1: polynomial, 2: radial, 3: sigmoid, 4: precomputed
        },
        "numfold" : 5,
        "corepos" : "center",
        "outdir":"outdir",
        "logit" : True
    }
