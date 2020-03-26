'''
1. Use all features, find best param
2. find the best top n
3. use those features  to find best param again
'''
import pandas as pd
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.model_selection import StratifiedKFold
import itertools

class BestModel:
    def __init__(self, clf, param_dict, train_data, topn = -1):
        self.classifier = clf
        self.param_dict = param_dict
        self.topn = topn
        self.train_data = train_data

        if self.topn == -1:
            self.topn = len(train_data.columns)-1

        self.cv_fold_num = 10
        self.cv_num = 1

    def run_all(self):
        # initialize classifier
        self.init_data(self.train_data)
        # get best hyperparam combination
        clf = self.get_best_param()
        new_x = self.train_data
        # get best top n
        if self.topn < len(self.train_data.columns) - 1:
            new_x = self.set_topn()
            # get best hyperparam for top n
            clf = self.get_best_param()
        # return the best model for top n
        return new_x, clf

    def init_data(self, df):
        self.x_train = df.loc[:, df.columns != 'label']
        self.y_train = df.loc[:, 'label']

    def set_clf(self, comb):
        if self.classifier == "RF":
            self.set_rf(comb)
        elif self.classifier == "DT":
            self.set_dt(comb)
        # d = {'RF': self.set_rf(comb),
        #      'DT': self.set_dt(comb)}
        # d[self.classifier]

    def set_rf(self, comb):
        self.clf = RandomForestClassifier(n_estimators=comb[0],
                                          max_depth=comb[1])

    def set_dt(self, comb):
        self.clf = DecisionTreeClassifier(criterion=comb[0],
                                        min_samples_split=comb[1],
                                        min_samples_leaf=comb[2])

    def set_lasso(self, comb):
        pass

    def set_mlp(self, comb):
        pass

    def set_svm(self, comb):
        pass

    def set_nb(self, comb):
        pass

    def get_best_param(self, save_to_file=False):
        # find the best model for the data
        x_train = self.x_train.values
        y_train = self.y_train.values


        comb_lst = list(self.param_dict.values())
        combinations = list(itertools.product(*comb_lst))
        keys  = list(self.param_dict.keys())
        output = []

        # loop through every hyperparameter combination
        max_auc = 0
        for comb in tqdm(combinations):
            self.set_clf(comb)
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
                best_comb = comb
                max_auc = avg_auc

        print("Best params:", best_comb)

        self.set_clf(best_comb)

        return self.clf

    def set_topn(self):
        model = self.clf.fit(self.x_train,self.y_train)
        feat_impt = model.feature_importances_
        print(sorted(zip(map(lambda x: round(x, 4), feat_impt), self.x_train.columns),
             reverse=True)[:self.topn])
        # get the new x
        self.x_train = self.train_data.iloc[:, feat_impt.argsort()[::-1][:self.topn]]
        df = self.x_train
        df['label'] = self.y_train
        return df

