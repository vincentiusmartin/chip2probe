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
    def __init__(self, clf, param_dict, train_data, topn = len(train_data.columns)-1):
        self.classifier = clf
        self.param_dict = param_dict
        self.topn = topn
        self.train_data = train_data

        self.cv_fold_num = 2
        self.cv_num = 1

    def run_all(self):
        # initialize classifier
        self.init_data(self.train_data)
        # get best hyperparam combination
        self.get_best_param()
        # get best top n 
        new_x = self.set_topn()
        # get best hyperparam for top n
        clf = self.get_best_param()
        # return the best model for top n
        return clf, new_x

    def init_data(self, df):
        self.x_train = df.loc[:, df.columns != 'label']
        self.y_train = df.loc[:, 'label']

    def set_clf(self, comb):
        d = {'RF': self.set_rf(comb),
             'DT': self.set_dt(comb)}
        d[self.classifier]

    def set_rf(self, comb):
        self.clf = RandomForestClassifier(n_estimators=comb[0],
                                          max_depth=comb[1])

    def set_dt(self, comb):
        self.clf = DecisionTreeClassifier()

    def set_lasso(self, comb):
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

        print("Best params:", best_comb[0], best_comb[1])

        self.set_clf(best_comb)

        return self.clf

    def set_topn(self):
        model = self.clf.fit(self.x_train,self.y_train)
        feat_impt = model.feature_importances_
        print(sorted(zip(map(lambda x: round(x, 4), feat_impt), self.x_train.columns), 
             reverse=True)[:self.topn])
        # get the new x 
        self.x_train = self.train_data.iloc[:, feat_impt.argsort()[::-1][:self.topn]]
        return self.x_train

    def get_topn(self):

        

if __name__=="__main__":


    # set hyperparameters to tune
    n_estimators = [i for i in range(2,21)]
    max_depth = [i for i in range(100,2001,100)]
    #n_estimators = [100]
    #max_depth = [10]
    comb_lst = [n_estimators, max_depth]
    

    # initialize an empty list for outputs
    output = []

    # set random state so the folds are the same for
    # each combination
    
    output = []
     
    #-------------------------Step 1------------------------
    
    '''
    #----------------------------Step 2--------------------------------

    # for the best set of parameters found above, find the top n features

    # get the x and y
    x_train = data.loc[:, data.columns != 'label']
    y_train = data.loc[:, 'label']
    output = []
    n_est = 20
    max_depth = 2000
    rf = RandomForestClassifier(n_estimators=n_est, max_depth=max_depth)
    model = rf.fit(x_train,y_train)
    # get the cv results for each of the top n features
    feat_impt = model.feature_importances_

    # for each top_n, what is the best param?
    # set random state so the folds are the same for
    # each combination
    cv = StratifiedKFold(n_splits=10, random_state = 19357, shuffle=True)
    for i in tqdm(range(1, len(feat_impt)+1)):
    #for i in tqdm(range(5, 6)):
        # get the new x and y values
        x = x_train.iloc[:, feat_impt.argsort()[::-1][:i]].values
        y = y_train.values
        max_auc = 0
        for comb in tqdm(combinations):
            n_est = comb[0]
            max_depth = comb[1]

            # initialize model
            rf = RandomForestClassifier(n_estimators=n_est, max_depth=max_depth)
            model = rf.fit(x, y)
            # perform 10 fold cross validation
            tot_acc = 0
            tot_auc = 0
            for train, test in cv.split(x, y):
                # fit the model
                model = rf.fit(x[train],y[train])
                predict_label = model.predict(x[test])
                predict_proba = model.predict_proba(x[test])[:,1]
                tot_acc += accuracy_score(y[test], predict_label)
                tot_auc += roc_auc_score(y[test], predict_proba)
            acc = tot_acc/10
            auc = tot_auc/10
            if auc > max_auc:
                best_comb = comb
                curr_acc, max_auc = acc, auc

        output.append([i, n_est, max_depth, curr_acc, max_auc])

        # write preliminary output
        columns = ["top_n", "best_n_est", "best_max_depth", "acc", "auc"] 
        output_df = pd.DataFrame.from_records(output, columns = columns)
        output_df.to_csv("rf_best_for_each_topn.csv", index=False)
    
    
    #-------------------------Step 3------------------------
    # using the top_n features found above, find the best params

    

    # loop through every hyperparameter combination
    max_auc = 0
    for comb in tqdm(combinations):
        n_est = comb[0]
        max_depth = comb[1]

        # initialize model
        rf = RandomForestClassifier(n_estimators=comb[0], max_depth=comb[1])

        # perform 10 fold cross validation
        tot_acc = 0
        tot_auc = 0
        for train, test in cv.split(x_train, y_train):
            # fit the model
            model = rf.fit(x_train[train],y_train[train])
            # get predictions
            predict_label = model.predict(x_train[test])
            predict_proba = model.predict_proba(x_train[test])[:,1]
            # calculate metrics
            tot_acc += accuracy_score(y_train[test], predict_label)
            tot_auc += roc_auc_score(y_train[test], predict_proba)
        # calculate average metrics
        acc = tot_acc/10
        auc = tot_auc/10
        output.append([comb[0], comb[1], acc, auc])
        # write preliminary output
        columns = ["n_est", "max_depth", "acc", "auc"] 
        output_df = pd.DataFrame.from_records(output, columns = columns)
        output_df.to_csv("rf_best_params_results_step3.csv", index=False)

        if auc > max_auc:
            curr_acc, max_auc = acc, auc
            best_comb = comb

    print("Best params:", best_comb[0], best_comb[1])
    '''
    print("DONE!!!")
    