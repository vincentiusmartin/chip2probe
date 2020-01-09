import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import itertools
import scipy

from sklearn import tree
from sklearn import metrics
from sklearn import ensemble
from sklearn import svm
from sklearn import linear_model
from sklearn import naive_bayes
from sklearn import model_selection
from sklearn import preprocessing
from matplotlib.backends.backend_pdf import PdfPages

from teacup.training import simpleclassifier

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
        if type(trainingdata) == str:
            self.training = pd.read_csv(trainingdata)
        elif type(trainingdata) == pd.core.frame.DataFrame:
            self.training = trainingdata[['sequence', 'bpos1', 'bpos2', 'distance', 'label']]
        self.motiflen = motiflen

    # ======== Getter ========

    def get_seq(self,tofile=False):
        seqdict = self.training['sequence'].to_dict()
        if not tofile:
            return seqdict
        else:
            keys = sorted(seqdict.keys())
            with open("sequences.txt",'w') as f:
                for key in keys:
                    f.write(">%s\n"%key)
                    f.write("%s\n"%seqdict[key])

    def get_labels_indexes(self):
        return self.training.groupby("label").groups

    # ======== Modifier to training data ========

    def get_numeric_label(self):
        train = self.training['label'].map({'cooperative': 1, 'additive': 0})
        return train

    # ======== Plot related ========

    def scatter_boxplot_col(self, colname, filepath="scatterbox.png"):
        groupdict = self.training.groupby(['label'])[colname].apply(list).to_dict()

        keys = groupdict.keys()
        listrep = [groupdict[key] for key in keys]

        pos = np.linspace(1,1+len(listrep)*0.5-0.5,len(listrep))
        bp = plt.boxplot(listrep,positions=pos,widths=0.4)
        plt.xticks(pos, keys)
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['caps'], color='black')

        for i in range(0,len(listrep)):
            y = listrep[i]
            x = np.random.normal(1+i*0.5, 0.02, size=len(y))
            plt.plot(x, y, 'r.', alpha=0.4,c='red')

        #print("Save distribution of row %s to %s" % (rownum,plotfilename))
        plt.savefig(filepath,positions=[0, 1])
        plt.clf() # clear canvas

    # ======== Training and testing modelss ========

    def extract_positional(self,seq):
        nucleotides = ['A','C','G','T']
        feature = []
        for k in range(1,3):
            perm = ["".join(p) for p in itertools.product(nucleotides, repeat=k)]
            for i in range(0,len(seq)+1-k):
                for kmer in perm:
                    if seq[i:i+k] == kmer:
                        feature.append(1)
                    else:
                        feature.append(0)
        return np.asarray(feature)

    def extract_kmer_ratio(self,seq):
        nucleotides = ['A','C','G','T']
        ratio = []
        for k in range(1,3):
            kmer = {}
            perm = ["".join(p) for p in itertools.product(nucleotides, repeat=k)]
            total = 0
            for p in perm:
                kmer[p] = 0
            for i in range(0,len(seq)+1-k):
                kmer[seq[i:i+k]] += 1
                total += 1
            for p in sorted(perm):
                ratio.append(float(kmer[p])/total)
        return ratio

    def extract_positional_features_bpos(self, seq, bpos1, bpos2, span_out=3, span_in=4):
        nucleotides = ['A','C','G','T']

        pos1 = bpos1 - 1
        pos2 = bpos2 - 1

        bseq1 = seq[pos1 - span_out : pos1 + span_in]
        bseq2 = seq[pos2 - span_in : pos2 + span_out]

        feature = self.extract_positional(bseq1) + self.extract_positional(bseq2)
        return feature

    def get_features(self,type="distance-numeric"):
        """
        type:
            distance-numeric
            distance-categorical
            sites-centered
            linker
        """
        if type == "distance-numeric":
            return self.training["distance"].values.reshape((-1,1))
        elif type == "distance-categorical":
            one_hot = pd.get_dummies(self.training['distance'])
            return  one_hot.values.tolist()
        elif type == "sites-centered":
            features = []
            for idx,row in self.training.iterrows():
                rowfeature = self.extract_kmer_features_bpos(row["sequence"],row["bpos1"],row["bpos2"])

                linker = row["sequence"][row["bpos1"] + self.motiflen // 2 : row["bpos2"] - self.motiflen // 2]
                ratio = self.extract_kmer_ratio(linker)

                all = np.concatenate((rowfeature,ratio,[self.training['distance'][idx]]))
                #features.append(preprocessing.normalize([all])[0])
                features.append(all)
            return features
        elif type == "sites-linker":
            features = []
            for idx,row in self.training.iterrows():
                numericdist = self.training["distance"].values.reshape((-1,1))
                # since the binding pos is one index, we need to -1
                midpos = (row["bpos2"] + row["bpos1"] - 1)//2
                seq = row["sequence"][midpos-13:midpos+13]
                features.append(self.extract_kmer_binary(seq) + [self.training['distance'][idx]])
            return features

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

    # ===========================

    def scatter_boxplot_dict(self, groupdict, filepath="scatterbox.png"):
        keys = groupdict.keys()
        listrep = [groupdict[key] for key in keys]

        pos = np.linspace(1,1+len(listrep)*0.5-0.5,len(listrep))
        bp = plt.boxplot(listrep,positions=pos,widths=0.4)
        plt.xticks(pos, keys)
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['caps'], color='black')

        for i in range(0,len(listrep)):
            y = listrep[i]
            x = np.random.normal(1+i*0.5, 0.02, size=len(y))
            plt.plot(x, y, 'r.', alpha=0.4,c='red')

        #print("Save distribution of row %s to %s" % (rownum,plotfilename))
        plt.savefig(filepath,positions=[0, 1])
        plt.clf() # clear canvas

    def boxplot_auc_distance(self):
        clfs = {
            #"decision tree":tree.DecisionTreeClassifier(),
            "random forest":ensemble.RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0),
            #"SVM":svm.SVC(kernel="rbf",gamma=1.0/5,probability=True),
            #"log regression":linear_model.LogisticRegression(),
            "simple":simpleclassifier.Simple1DClassifier(),
            #"gradient boosting":ensemble.GradientBoostingClassifier(),
            #"naive bayes":naive_bayes.GaussianNB()
           }

        dists = ["distance-numeric","distance-categorical"]

        auc_dict = {}
        for dist in dists:
            auc_dict[dist] = []
            for i in range(100):
                x_train = self.get_features(dist)
                y_train = self.get_numeric_label().values
                fpr_list, tpr_list, auc_list = self.test_with_cv(clfs, x_train, y_train)
                auc_dict[dist].append(auc_list[0])

        self.scatter_boxplot_dict(auc_dict)


    def get_features_custom(self, span_out=3, span_in=4, disttype="distance-categorical"):
        features = []
        dist = ""
        if disttype == "distance-categorical":
            dist = pd.get_dummies(self.training['distance']).values.tolist()
        elif disttype == "distance-numeric":
            dist = self.training['distance']

        for idx,row in self.training.iterrows():
            rowfeature = self.extract_positional_features_bpos(row["sequence"],row["bpos1"],row["bpos2"], span_out, span_in)

            linker = row["sequence"][row["bpos1"] + self.motiflen // 2 : row["bpos2"] - self.motiflen // 2]
            ratio = self.extract_kmer_ratio(linker)

            all = np.concatenate((rowfeature,ratio))
            if disttype == "distance-categorical":
                all = np.concatenate((all,dist[idx]))
            elif disttype == "distance-numeric":
                all = np.concatenate((all,[dist[idx]]))
            #features.append(preprocessing.normalize([all])[0])
            features.append(all)
        return features

    def test_seq_features(self,outpath="auc.png"):
        clfs = {
                "decision tree":tree.DecisionTreeClassifier(),
                "random forest":ensemble.RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0),
                #"SVM":svm.SVC(kernel="rbf",gamma=1.0/5,probability=True),
                #"log regression":linear_model.LogisticRegression(),
                "simple":simpleclassifier.Simple1DClassifier(),
                #"gradient boosting":ensemble.GradientBoostingClassifier(),
                #"naive bayes":naive_bayes.GaussianNB()
               }

        span_in_list = [1,2,3,4,5,6,7,8]
        span_out_list = [1,2,3]
        dlist = ["distance-categorical","distance-numeric","None"]
        combs = [{"distance":"distance-numeric"},{"distance":"distance-categorical"}]

        for x1 in span_in_list:
            for x2 in span_out_list:
                for x3 in dlist:
                    combs.append({"span_in":x1,"span_out":x2,"wdist":x3})

        classifier_names = list(clfs.keys())

        # we only need to make y_train once
        y_train = self.get_numeric_label().values

        n = 0
        numcol = 2 # adjust the number of columns in the plot
        numrow = 2 # adjust the number of rows in the plot
        with PdfPages(outpath) as pdf:
            for comb in combs:
                # we need this because n is always reset after one full page
                if n == 0:
                    fig = plt.figure(figsize=(12,12))
                    fig.subplots_adjust(hspace=0.4,wspace=0.5)
                n+=1

                if "distance" in comb:
                    x_train = self.get_features(comb["distance"])
                else:
                    x_train = self.get_features_custom(comb["span_out"],comb["span_in"],comb["wdist"])
                #self.display_output(fpr_list, tpr_list, auc_list, list(clfs.keys()), path=outpath)

                fpr_list, tpr_list, auc_list = self.test_with_cv(clfs, x_train, y_train)

                ax = fig.add_subplot(numcol,numrow,n)

                ax.plot([0, 1], [0, 1], linestyle="--", color="red", alpha=0.1)
                for i in range(len(fpr_list)):
                    ax.plot(fpr_list[i], tpr_list[i], lw=2, alpha=0.4, label='%s, AUC %f' % (classifier_names[i], auc_list[i]))

                # Show the ROC curves for all classifiers on the same plot
                ax.xaxis.set_label_text('False Positive Rate')
                ax.yaxis.set_label_text('True Positive Rate')
                if "distance" in comb:
                    ax.set_title("distance type %s" % comb["distance"])
                else:
                    ax.set_title('span_out %d, span_in %d, with_dist %s' % (comb["span_out"],comb["span_in"],comb["wdist"]))
                ax.legend(loc="lower right")
                if n == numcol*numrow:
                    pdf.savefig(fig)
                    plt.close()
                    n = 0
            pdf.savefig(fig)
            plt.close()

    # ===========================

    def test_model(self, feature_type, testing_type="cv", outpath="roc.png"):
        """
        testing_type:
            cv: cross validation
            train: test on train
        """

        x_train = self.get_features(feature_type)
        y_train = self.get_numeric_label().values
        #print(len(x_train),len(y_train))

        clfs = {
                "decision tree":tree.DecisionTreeClassifier(),
                "random forest":ensemble.RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0),
                #"SVM":svm.SVC(kernel="rbf",gamma=1.0/5,probability=True),
                #"log regression":linear_model.LogisticRegression(),
                #"simple":simpleclassifier.Simple1DClassifier(),
                #"gradient boosting":ensemble.GradientBoostingClassifier(),
                #"naive bayes":naive_bayes.GaussianNB()
               }

        if testing_type == "cv":
            fpr_list, tpr_list, auc_list = self.test_with_cv(clfs, x_train, y_train)
        else:
            fpr_list, tpr_list, auc_list = self.test_on_train(clfs,x_train,y_train)

        self.display_output(fpr_list, tpr_list, auc_list, list(clfs.keys()), path=outpath)

    def test_with_cv(self,clfs,x_train,y_train,fold=10):
        fpr_list = []
        tpr_list = []
        auc_list = []
         # Compute ROC curve and ROC area with averaging for each classifier
        for key in clfs:
            base_fpr = np.linspace(0, 1, 101)
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
                    fpr, tpr, _ = metrics.roc_curve(lbl_test, y_score[:, 1])
                    auc = metrics.roc_auc_score(lbl_test, y_score[:,1])
                    #print("fold " + str(i) + " AUC: " + str(auc))
                    # vmartin: please have the package name instead of using
                    # the function directly so we know where does the function
                    # come from :)
                    tpr = scipy.interp(base_fpr, fpr, tpr)
                    tprs.append(tpr)
                    aucs_val.append(auc)
                    i += 1

            # calculate mean true positive rate
            tprs = np.array(tprs)
            mean_tprs = tprs.mean(axis=0)

            # calculate mean auc
            aucs_val = np.array(aucs_val)
            mean_aucs = aucs_val.mean(axis=0)

            fpr_list.append(base_fpr)
            tpr_list.append(mean_tprs)
            auc_list.append(mean_aucs)

        return fpr_list, tpr_list, auc_list

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


    def display_output(self, fpr_list, tpr_list, auc_list, classifier_names, path):
        """
            This plots the average ROC curve of all the classifiers in a single plot
        """
        plt.clf() # first, clear the canvas

        plt.plot([0, 1], [0, 1], linestyle="--", color="red", alpha=0.1)
        for i in range(len(fpr_list)):
            plt.plot(fpr_list[i], tpr_list[i], lw=2, alpha=0.4, label='%s, AUC %f' % (classifier_names[i], auc_list[i]))

        # Show the ROC curves for all classifiers on the same plot
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Average ROC Curves for All Classifiers')
        plt.legend(loc="lower right")
        plt.savefig(path)
