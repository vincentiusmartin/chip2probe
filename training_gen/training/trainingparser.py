import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import itertools
import scipy
import subprocess

from sklearn import tree
from sklearn import metrics
from sklearn import ensemble
from sklearn import svm
from sklearn import linear_model
from sklearn import naive_bayes
from sklearn import model_selection
from sklearn import preprocessing

from teacup.training import simpleclassifier
from teacup.training import dnashape
from teacup import utils

# TODO: hard coded for now, make it better later
shapepath = "/Users/vincentiusmartin/Research/Cooperativity/teacup/data/dnashape"

def calculate_fpr_tpr(ytrue,ypred):
    if len(ytrue) != len(ypred):
        print("the length of y-true and y-pred differ")
        return 0
    fp_count = 0
    tp_count = 0
    pos_count = 0
    neg_count = 0
    for i in range(len(ytrue)):
        if ytrue[i] == 1:
            pos_count += 1
            if ypred[i] == 1:
                tp_count += 1
        elif ytrue[i] == 0:
            neg_count += 1
            if ypred[i] == 1:
                fp_count += 1
    fpr = float(fp_count)/neg_count
    tpr = float(tp_count)/pos_count
    return fpr,tpr

class TrainingParser:

    def __init__(self, trainingdata,motiflen):
        if type(trainingdata) == str: # input path
            self.training = pd.read_csv(trainingdata)
        elif type(trainingdata) == pd.core.frame.DataFrame: # from an existing data frame
            self.training = trainingdata[['sequence', 'bpos1', 'bpos2', 'distance', 'label']]
        self.motiflen = motiflen

    # ===== Getter part ====

    def get_labels_indexes(self):
        return self.training.groupby("label").groups

    def get_bsites(self):
        '''
        return a list of bpos1,bpos2 where each is a dictionary
        '''
        sitedict = {}
        bpos1 = self.training["bpos1"].to_dict()
        bpos2 = self.training["bpos2"].to_dict()
        return bpos1,bpos2

    def get_seq(self,tofile=False):
        seqdict = self.training["sequence"].to_dict()
        if not tofile:
            return seqdict
        else:
            keys = sorted(seqdict.keys())
            with open("sequences.txt",'w') as f:
                for key in keys:
                    key_int = int(key)
                    f.write(">%d\n"%key_int)
                    f.write("%s\n"%seqdict[key])

    def get_seq_aligned(self,tofile=False):
        '''
        this function align sequence based on its first bpos location, useful
        to get shape features
        '''
        min_dist = self.training.min()["bpos1"]
        max_dist = self.training.max()["bpos1"]
        trimlen = len(self.training.iloc[0]["sequence"]) - (max_dist - min_dist)

        align_dict = {}
        pr = True
        for idx,row in self.training.iterrows():
            lshift = row["bpos1"] - min_dist
            if pr:
                print("Binding pos location, b1: %d, b2: %d" % (min_dist,row["bpos2"]-lshift))
                pr = False
            aligned = row["sequence"][lshift:]
            trimmed = row["sequence"][:trimlen]
            align_dict[idx] = trimmed
        if not tofile:
            return align_dict
        else:
            return utils.dictlist2file(align_dict,"sequence.txt")

    # ========= test model =========
    def test_model(self, feature_list, testing_type="cv", outpath="roc.png"):
        """
        testing_type:
            cv: cross validation
            train: test on train
            Produce AUC
        """

        x_train = self.get_features(feature_list)
        y_train = self.get_numeric_label().values
        #print(len(x_train),len(y_train))

        clfs = {
                #"decision tree":tree.DecisionTreeClassifier(),
                "random forest":ensemble.RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0),
                #"SVM":svm.SVC(kernel="rbf",gamma=1.0/5,probability=True),
                #"log regression":linear_model.LogisticRegression(),
                "simple":simpleclassifier.Simple1DClassifier(),
                #"gradient boosting":ensemble.GradientBoostingClassifier(),
                #"naive bayes":naive_bayes.GaussianNB()
               }

        if testing_type == "cv":
            fpr_list, tpr_list, auc_list = self.test_with_cv(clfs, x_train, y_train)
        else:
            fpr_list, tpr_list, auc_list = self.test_on_train(clfs,x_train,y_train)

        self.display_output(fpr_list, tpr_list, auc_list, path=outpath)

# ========= Visualization =========

    def visualize_random_forest(self, types, max_depth=10):
        rf = ensemble.RandomForestClassifier(n_estimators=100, max_depth=max_depth,random_state=0)
        # print trees
        # feature importance
        x_df = self.get_features(types,ret_tbl=True)
        x_train = x_df.values.tolist()

        y_train = self.get_numeric_label().values
        rf.fit(x_train,y_train)

        # draw a tree from the forest, let's say tree 5
        estimator = rf.estimators_[5]
        tree.export_graphviz(estimator, out_file='tree.dot',
                feature_names = x_df.columns,
                class_names = ['additive','cooperative'],
                rounded = True, proportion = False,
                precision = 2, filled = True)
        subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])

        # do feature importance, code is taken from Farica's
        feature_importances = pd.DataFrame(rf.feature_importances_,
                                           index = x_df.columns,
                                           columns=['importance']).sort_values('importance',ascending=False)
        print(feature_importances)


    def get_features(self, types, ret_tbl=False):
        """
        type:
            dist-numeric
            dist-categorical
            linker_[k]mer
            positional_in_[x]_out_[y]
            shape
        ret_tbl:
            False: return a list of list
            True: return a list of dictionary--this can be directly converted
                  into a data frame.
        """
        if not (isinstance(types, list) or isinstance(types, tuple)):
            print("Error: Input types must be a list or a tuple!")
            return []

        features = []
        for feature_type in types:
            if feature_type == "dist-numeric":
                # (self.training["distance"].values.reshape((-1,1)))
                rfeature = [{"dist-numeric":x} for x in self.training["distance"].values]
            elif feature_type == "dist-categorical":
                one_hot = pd.get_dummies(self.training['distance'])
                one_hot.columns = ["dist-num-%d"%col for col in one_hot.columns]
                #features.append(one_hot.values.tolist())
                rfeature = one_hot.to_dict('records')
            elif feature_type.startswith("linker"):
                # this uses kmer ratio
                rfeature = []
                for idx,row in self.training.iterrows():
                    start = row["bpos1"] + self.motiflen // 2 + 1
                    end = row["bpos2"] - self.motiflen // 2
                    linker = row["sequence"][start:end]
                    k = int(feature_type[len("linker_") : feature_type.find("mer")])
                    ratio = utils.extract_kmer_ratio(linker,k)
                    rfeature.append(ratio)
                    #ratio_feature = [x[1] for x in sorted(ratio.items(), key=lambda k:k[0])]
                    #rowfeatures.append(ratio_feature)
            elif feature_type.startswith("positional"):
                splitted = feature_type.split("_")
                s_in = int(splitted[2])
                s_out = int(splitted[4])
                rfeature = []
                for idx,row in self.training.iterrows():
                    pos_feature = utils.extract_positional_features(row["sequence"], row["bpos1"], row["bpos2"],
                                 span_out=s_out, span_in=s_in)
                    rfeature.append(pos_feature)
            elif feature_type == "shape":
                ds = dnashape.DNAShapes(shapepath,self.get_bsites())
                rfeature = ds.get_features()

            features = utils.merge_listdict(features,rfeature)

        df_features = pd.DataFrame(features)
        if ret_tbl: # return as data frame
            return df_features
        else:
            return df_features.values.tolist()

    # ======== Modifier to training data ========

    def get_numeric_label(self):
        train = self.training['label'].map({'cooperative': 1, 'additive': 0})
        return train

    # ======= For simple model that is based on distance only =======
    def roc_simple_clf(self,n_splits=1):
        # still numeric for now
        x_train = self.training["distance"].values
        y_train = self.get_numeric_label().values
        distances = self.training['distance'].unique()

        if n_splits > 1:
            cv = model_selection.KFold(n_splits=n_splits,shuffle=True)
            split = cv.split(x_train,y_train)
        else:
            split = [(range(len(x_train)),range(len(y_train)))]

        fpr_all = []
        tpr_all = []
        auc_all = []

        for train, test in split:
            fpr_list = [0]
            tpr_list = [0]
            for dist in sorted(distances):
                scf = simpleclassifier.Simple1DClassifier()
                scf.fit_on_thres(x_train[train],y_train[train],dist)
                y_pred = scf.test(x_train[test])
                #print("Accuracy %f" % metrics.accuracy_score(ytrain, ypred))
                fpr,tpr = calculate_fpr_tpr(y_train[test], y_pred)
                fpr_list.append(fpr)
                tpr_list.append(tpr)

            fpr_list.append(1)
            tpr_list.append(1)

            auc = metrics.auc(fpr_list,tpr_list)
            auc_all.append(auc)
            fpr_all.append(fpr_list)
            tpr_all.append(tpr_list)
        return fpr_all,tpr_all,auc_all

    # ====== Processing part ======

    def compare_distance_features(self, iter=10, fpr_lim=100):
        clfs = {
            #"decision tree":tree.DecisionTreeClassifier(),
            "random forest":ensemble.RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0),
            #"SVM":svm.SVC(kernel="rbf",gamma=1.0/5,probability=True),
            #"log regression":linear_model.LogisticRegression(),
            "simple":simpleclassifier.Simple1DClassifier(),
            #"gradient boosting":ensemble.GradientBoostingClassifier(),
            #"naive bayes":naive_bayes.GaussianNB()
           }

        dists = [["dist-numeric"],["dist-categorical"]]

        auc_dict = {}
        for dist_type in dists:
            dname = dist_type[0]
            auc_dict[dname] = []
            for i in range(iter):
                print("Processing using %s, iteration %d" % (dist_type,i+1))
                x_train = self.get_features(dist_type)
                y_train = self.get_numeric_label().values
                fpr_list, tpr_list, auc_list = self.test_with_cv(clfs, x_train, y_train,fpr_lim=fpr_lim)
                auc_dict[dname].append(auc_list['random forest'])


        print("Making scatter boxplot for each feature...")
        utils.scatter_boxplot_dict(auc_dict,ylabel="AUC")

        print("Two sided wilcox test, pval: %.4f" % utils.wilcox_test(auc_dict["dist-numeric"],auc_dict["dist-categorical"]))
        print("Numeric > Categorical test, pval: %.4f" % utils.wilcox_test(auc_dict["dist-numeric"],auc_dict["dist-categorical"],alternative="greater"))
        print("Numeric < Categorical test, pval: %.4f" % utils.wilcox_test(auc_dict["dist-numeric"],auc_dict["dist-categorical"],alternative="less"))

    def compare_dist_pos_features(self, iter=10, fpr_lim=100, path="dist_positional.pdf"):
        clfs = {
            "random forest":ensemble.RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0)
        }
        y_train = self.get_numeric_label().values

        span_out_list = [0,1,2,3]
        span_in_list = [0,1,2,3,4,5,6,7,8]
        #spans = list(itertools.product(span_in_list, span_out_list))

        auc_all = []
        for so in span_out_list:
            auc_dict = {}
            for si in span_in_list:
                type = "positional_in_%d_out_%d" % (si,so)
                #print(type)
                features = [type,"dist-numeric"]
                x_train = self.get_features(features)
                fea_name = ",".join(features)
                auc_dict[fea_name] = []
                for i in range(iter):
                    fpr_list, tpr_list, auc_list = self.test_with_cv(clfs, x_train, y_train,fpr_lim=fpr_lim)
                    auc_dict[fea_name].append(auc_list['random forest'])
            auc_all.append(auc_dict)
            #x_train = self.get_features([])
        utils.multiple_scatter_boxplots(auc_all,ylabel="AUC",filepath=path)


    # super important function
    def compare_prefix_features(self, features, iter=10, fpr_lim=100, max_depth=10, max_comb = 2, path="linker.png"):
        #prefix = ["dist-numeric", "linker_1mer", "linker_2mer"]
        y_train = self.get_numeric_label().values

        clfs = {
            "random forest":ensemble.RandomForestClassifier(n_estimators=100, max_depth=max_depth)
        }

        auc_dict = {}
        for i in range(max_comb):
            for comb in itertools.combinations(features, i+1):
                comb_name = ", ".join(comb)
                auc_dict[comb_name] = []
                for i in range(iter):
                    print("Processing using %s, iteration %d" % (str(comb_name),i+1))
                    x_train = self.get_features(comb)
                    fpr_list, tpr_list, auc_list = self.test_with_cv(clfs, x_train, y_train,fpr_lim=fpr_lim)
                    auc_dict[comb_name].append(auc_list['random forest'])
        utils.scatter_boxplot_dict(auc_dict,ylabel="AUC",filepath=path)

        keys = auc_dict.keys()
        for comb in itertools.combinations(keys, 2):
            print("Two sided wilcox test, pval: %.4f" % utils.wilcox_test(auc_dict[comb[0]], auc_dict[comb[1]]))
            print("%s > %s test, pval: %.4f" % (comb[0],comb[1],utils.wilcox_test(auc_dict[comb[0]], auc_dict[comb[1]], alternative="greater")) )
            print("%s > %s test, pval: %.4f" % (comb[1],comb[0],utils.wilcox_test(auc_dict[comb[0]], auc_dict[comb[1]], alternative="less")) )
            print("---------------------")

    def test_with_cv(self,clfs,x_train,y_train,fold=10,fpr_lim=100):
        fpr_dict = {}
        tpr_dict = {}
        auc_dict = {}
         # Compute ROC curve and ROC area with averaging for each classifier
        for key in clfs:
            # we limit this to get roc curve / auc until the fpr that we want
            base_fpr = np.linspace(0, 1, 101)[:fpr_lim+1]
            tprs = []
            aucs_val = []
            if key == "simple":
                fprs_simple,tprs_simple,aucs_val = self.roc_simple_clf(n_splits=fold)
                for i in range(0,len(fprs_simple)):
                    tpr = scipy.interp(base_fpr, fprs_simple[i], tprs_simple[i])
                    tprs.append(tpr)
            else:
                cv = model_selection.KFold(n_splits=fold,shuffle=True)
                # initialize a list to store the average fpr, tpr, and auc
                print("Cross validation on %s" % key)
                i = 1
                for train_idx,test_idx in cv.split(x_train,y_train):
                    # need to convert this with index, somehow cannot do
                    # x_train[train_idx] for multi features
                    data_train = [x_train[i] for i in train_idx]
                    data_test = [x_train[i] for i in test_idx]
                    lbl_train = [y_train[i] for i in train_idx]
                    lbl_test = [y_train[i] for i in test_idx]

                    model = clfs[key].fit(data_train, lbl_train)
                    y_score = model.predict_proba(data_test)
                    #print(len(y_score),len(train_idx),len(test_idx))
                    print(lbl_test,y_score[:, 1])
                    fpr, tpr, _ = metrics.roc_curve(lbl_test, y_score[:, 1])
                    #auc = metrics.roc_auc_score(lbl_test, y_score[:,1])
                    tpr = scipy.interp(base_fpr, fpr, tpr) # add points to the plotting
                    res_auc = metrics.auc(base_fpr, tpr)
                    tprs.append(tpr)
                    aucs_val.append(res_auc)
                    i += 1
                    #break

            # calculate mean true positive rate
            tprs = np.array(tprs)
            mean_tprs = tprs.mean(axis=0)

            # calculate mean auc
            aucs_val = np.array(aucs_val)
            mean_aucs = aucs_val.mean(axis=0)

            fpr_dict[key] = base_fpr
            tpr_dict[key] = mean_tprs
            auc_dict[key] = mean_aucs

        return fpr_dict, tpr_dict, auc_dict

    def test_on_train(self,clfs,x_train,y_train):
        auc_total = 0
        fpr_list = []
        tpr_list = []
        auc_list = []
        for key in clfs:
            if key == "simple":
                fpr,tpr,auc = self.roc_simple_clf()
                fpr = fpr[0]
                tpr = tpr[0]
                auc = auc[0]
                #plt.plot(fpr,tpr,label="distance threshold, training auc=%f" % auc,linestyle=":", color="orange")
            else:
                print("key is:", key)
                clf = clfs[key].fit(x_train, y_train)
                y_pred = clf.predict_proba(x_train)[:, 1]

                # https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
                # print("Accuracy %s: %f" % (key,metrics.accuracy_score(y_train, y_pred)))

                # ROC curve
                fpr, tpr, _ = metrics.roc_curve(y_train, y_pred)
                auc = metrics.roc_auc_score(y_train, y_pred)
                #plt.plot(fpr,tpr,label="%s, training auc=%f" % (key,auc))

            fpr_list.append(fpr)
            tpr_list.append(tpr)
            auc_list.append(auc)
            auc_total += auc
        print("Average AUC %f"%(auc_total/len(clfs)))

        return fpr_list, tpr_list, auc_list


    # ========= Plotting =========

    def display_output(self, fpr_dict, tpr_dict, auc_dict, path):
        """
            This plots the average ROC curve of all the classifiers in a single plot
        """
        plt.clf() # first, clear the canvas

        plt.plot([0, 1], [0, 1], linestyle="--", color="red", alpha=0.1)
        for key in fpr_dict:
            plt.plot(fpr_dict[key], tpr_dict[key], lw=2, alpha=0.4, label='%s, AUC %f' % (key, auc_dict[key]))

        # Show the ROC curves for all classifiers on the same plot
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Average ROC Curves for All Classifiers')
        plt.legend(loc="lower right")
        plt.savefig(path)
