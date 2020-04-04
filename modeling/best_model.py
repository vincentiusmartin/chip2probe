'''
1. Use all features, find best param
2. find the best top n
3. use those features  to find best param again
'''
import pandas as pd
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
import sklearn.metrics as metrics
from sklearn.model_selection import StratifiedKFold
import itertools
import concurrent.futures as cc
import timeit
import os
import functools

class BestModel:
    def __init__(self, clf, param_dict, train_data, topn = -1, cv_fold_num=10, cv_num=1):
        self.classifier = clf
        self.param_dict = param_dict
        self.topn = topn # integer
        self.train_data = train_data

        if self.topn == -1:
            self.topn = len(train_data.columns)-1

        self.cv_fold_num = cv_fold_num
        self.cv_num = cv_num

    def run_all(self, num_workers=os.cpu_count(), score_type="auc"):
        # initialize classifier
        self.init_data(self.train_data)
        # get best hyperparam combination
        clf = self.get_best_param(num_workers=num_workers,score=score_type)
        new_x = self.train_data
        # get best top n
        if self.topn < len(self.train_data.columns) - 1:
            new_x = self.set_topn(clf)
            # get best hyperparam for top n
            clf = self.get_best_param(num_workers=num_workers,score=score_type)
        # return the best model for top n
        return new_x, clf

    def init_data(self, df):
        self.x_train = df.loc[:, df.columns != 'label']
        self.y_train = df.loc[:, 'label']

    def get_clf(self, comb_dict):
        if self.classifier == "RF":
            return RandomForestClassifier(**comb_dict)
        elif self.classifier == "DT":
            return DecisionTreeClassifier(**comb_dict)
        # d = {'RF': self.set_rf(comb),
        #      'DT': self.set_dt(comb)}
        # d[self.classifier]


    # def set_lasso(self, comb):
    #     pass
    #
    # def set_mlp(self, comb):
    #     pass
    #
    # def set_svm(self, comb):
    #     pass
    #
    # def set_nb(self, comb):
    #     pass

    def run_kfold(self, comb, score="auc"):
        """
        metrics: auc (Area Under ROC Curve)/pr (Precision Recall)
        """
        x_train = self.x_train.values
        y_train = self.y_train.values
        keys  = list(self.param_dict.keys())

        comb_dict = {keys[i]:comb[i] for i in range(len(keys))}
        clf = self.get_clf(comb_dict)
        avg_acc = 0
        avg_scr = 0
        # take the average of runs
        for i in range(self.cv_num):
            # perform 10 fold cross validation
            tot_acc = 0
            tot_score = 0
            cv = StratifiedKFold(n_splits=self.cv_fold_num, shuffle=True)
            for train, test in cv.split(x_train, y_train):
                # fit the model
                model = clf.fit(x_train[train],y_train[train])
                # get predictions
                predict_label = model.predict(x_train[test])
                predict_proba = model.predict_proba(x_train[test])[:,1]
                # calculate metrics
                tot_acc += metrics.accuracy_score(y_train[test], predict_label)
                if score == "pr":
                    tot_score += metrics.average_precision_score(y_train[test], predict_proba)
                elif score == "auc" or score == "acc":
                    # if accuracy we still return auc
                    tot_score += metrics.roc_auc_score(y_train[test], predict_proba)
                else:
                    raise Exception("score metrics should be pr/auc/acc")
            # calculate average metrics
            acc = tot_acc/self.cv_fold_num
            scr = tot_score/self.cv_fold_num
            avg_acc += acc
            avg_scr += scr
        avg_acc = avg_acc/self.cv_num
        avg_scr = avg_scr/self.cv_num
        return {"params":comb_dict, "avg_acc": avg_acc, "avg_%s"%score: avg_scr}

    def get_best_param(self, save_to_file=False, num_workers=os.cpu_count(), score="auc"):
        print("==== Get best param, utilizing %d workers =====" % num_workers)
        start_time = timeit.default_timer()
        comb_lst = list(self.param_dict.values())
        combinations = list(itertools.product(*comb_lst))

        # loop through every hyperparameter combination
        max_auc = 0
<<<<<<< HEAD
        for comb in tqdm(combinations):
            comb_dict = {keys[i]:comb[i] for i in range(len(keys))}
            self.set_clf(comb_dict)
            avg_acc = 0
            avg_auc = 0
            # take the average of runs
            for i in range(self.cv_num):
                # perform 10 fold cross validation
                tot_acc = 0
                tot_auc = 0
                cv = StratifiedKFold(n_splits=self.cv_fold_num, shuffle=True)
                for train, test in cv.split(x_train, y_train):
                    # fit the model
                    model = self.clf.fit(x_train[train],y_train[train])
                    # get predictions
                    predict_label = model.predict(x_train[test])
                    predict_proba = model.predict_proba(x_train[test])[:,1]
                    # calculate metrics
                    tot_acc += accuracy_score(y_train[test], predict_label)
                    tot_auc += roc_auc_score(y_train[test], predict_proba)
                # calculate average metrics
                acc = tot_acc/self.cv_fold_num
                auc = tot_auc/self.cv_fold_num
                avg_acc += acc
                avg_auc += auc
            avg_acc = avg_acc/self.cv_num
            avg_auc = avg_auc/self.cv_num
            output.append([comb[0], comb[1], avg_acc, avg_auc])
            if save_to_file:
                # write preliminary output
                columns = ["n_est", "max_depth", "acc", "auc"]
                output_df = pd.DataFrame.from_records(output, columns = columns)
                output_df.to_csv("rf_all_params_results_step1.csv", index=False)

            if avg_auc > max_auc:
                best_comb = comb_dict
                max_auc = avg_auc

        print("Best params:", best_comb)
        self.set_clf(best_comb)

        return self.clf

    def set_topn(self):
        model = self.clf.fit(self.x_train,self.y_train)
=======
        run_kfold_partial = functools.partial(self.run_kfold, score=score)
        with cc.ProcessPoolExecutor(max_workers = num_workers) as executor:
            # TODO: update input to combinations to dictionary
            output = list(tqdm(executor.map(run_kfold_partial, combinations), total=len(combinations)))

            # if save_to_file:
            #     # write preliminary output
            #     columns = ["n_est", "max_depth", "acc", "auc"]
            #     output_df = pd.DataFrame.from_records(output, columns = columns)
            #     output_df.to_csv("rf_all_params_results_step1.csv", index=False)

        best = max(output, key=lambda x:x['avg_%s'%score])

        print("Best params (using %s): %s" % (score,str(best["params"])))
        print("Total running time: %d seconds" % (timeit.default_timer() - start_time))
        return self.get_clf(best["params"])

    def set_topn(self, clf):
        model = clf.fit(self.x_train,self.y_train)
>>>>>>> 7ffdc5ef0f76ce847cdd021334dd7e197edd1876
        feat_impt = model.feature_importances_
        print("Top %d features"%self.topn,sorted(zip(map(lambda x: round(x, 4), feat_impt), self.x_train.columns),
             reverse=True)[:self.topn])
        # get the new x
        self.x_train = self.train_data.iloc[:, feat_impt.argsort()[::-1][:self.topn]]
        df = self.x_train
        df['label'] = self.y_train
        return df
