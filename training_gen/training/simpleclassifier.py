from sklearn import metrics


class Simple1DClassifier:
    """
    input table needs to have column label on it
    """

    def __init__(self):
        self.label_gt = -1 # greater than
        self.label_lt = -1 # less than
        self.threshold = 0

    def update_params(self,label_gt,label_lt,threshold):
        self.label_gt = label_gt
        self.label_lt = label_lt
        self.threshold = threshold

    def label_from_dist(self,x_train,y_train,threshold):
        index = [y_train[i] for i in range(len(x_train)) if x_train[i] >= threshold]
        label_gt = max(index,key=index.count)
        if label_gt == 1:
            label_lt = 0
        else:
            label_lt = 1
        threshold = threshold
        return label_gt,label_lt,threshold

    def predict_with_thres(self,xtest,label_gt,label_lt,threshold):
        predictions = []
        for x in xtest:
            if x >= threshold:
                predictions.append(label_gt)
            else:
                predictions.append(label_lt)
        return predictions

    def fit_on_thres(self,x_train,y_train,threshold):
        label_gt, label_lt, threshold = self.label_from_dist(x_train,y_train,threshold)
        self.update_params(label_gt, label_lt, threshold)

    def fit_best_thres(self, x_train, y_train, threslist):
        best_thres = -1
        best_acc = -1
        best_lgt = -1
        best_llt = -1
        for thres in threslist:
            label_gt,label_lt,threshold = self.label_from_dist(x_train,y_train,thres)
            #
            y_pred = self.predict_thres(x_train, label_gt, label_lt, threshold)
            acc = metrics.accuracy_score(y_train, y_pred)
            if acc > best_acc:
                best_acc = acc
                best_thres = thres
                best_lgt = label_gt
                best_llt = label_lt
        self.label_gt = best_lgt
        self.label_lt = best_llt
        self.threshold = best_thres

    def test(self,xtest):
        return self.predict_with_thres(xtest, self.label_gt, self.label_lt, self.threshold)
