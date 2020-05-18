# Modeling

## Generating cooperative model

### Useful import


### Reading input training data

The `flip_one_face_orientation` makes sure all sequences with HT/TH orientations
face the same way.
~~~~{.python}
# trainingpath = "input/training_data/training_p01_adjusted.tsv"
# df = pd.read_csv(trainingpath, sep="\t")
# traindf = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])
~~~~

### Model with dist, orientation, and imads score

~~~~{.python}
# dist_ori_imads = t.get_feature_all({
#     "distance":{"type":"numerical"},
#     "sitepref": {"imadsmodel": imads12, "modelwidth":12},
#     "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
# })
~~~~
