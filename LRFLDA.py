import os
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import scale,StandardScaler
from sklearn.metrics import roc_curve, auc
import utils.tools as utils
from sklearn.model_selection import LeaveOneOut
data_train1=pd.read_excel('TrainingSample_Lasso.xlsx',header = None)
data_=np.array(data_train1)
data_train11=np.array(data_train1)
data=data_[:,1:]
label=data_[:,0]
def get_shuffle(data,label):
    index = [i for i in range(len(label))]
    np.random.shuffle(index)
    data = data[index]
    label = label[index]
    return data,label
X_=scale(data)
X,y=get_shuffle(X_,label)
sepscores = []
y_score=np.ones((1,2))*0.5
y_class=np.ones((1,1))*0.5
loo = LeaveOneOut()
for train, test in loo.split(X):
       cv_clf = RandomForestClassifier(n_estimators=500, criterion='gini', max_depth=10,
                                    min_samples_split=2, min_samples_leaf=1,
                                    min_weight_fraction_leaf=0.0, max_features='auto',
                                    max_leaf_nodes=None, min_impurity_decrease=0.0,
                                    bootstrap=True,
                                    oob_score=False, n_jobs=1, random_state=None, verbose=0,
                                    warm_start=False, class_weight=None)
    X_train=X[train]
    y_train=y[train]
    X_test=X[test]
    y_test=y[test]
    y_sparse=utils.to_categorical(y)
    y_train_sparse=utils.to_categorical(y_train)
    y_test_sparse=utils.to_categorical(y_test)
    hist=cv_clf.fit(X_train, y_train)
    y_predict_score=cv_clf.predict_proba(X_test)
    y_predict_class= utils.categorical_probas_to_classes(y_predict_score)
    y_score=np.vstack((y_score,y_predict_score))
    y_class=np.vstack((y_class,y_predict_class))
    cv_clf=[]
y_class=y_class[1:]
y_score=y_score[1:]
fpr, tpr, _ = roc_curve(y_sparse[:,0], y_score[:,0])
roc_auc = auc(fpr, tpr)
acc, precision,npv, sensitivity, specificity, mcc,f1 = utils.calculate_performace(len(y_class), y_class, y)
result=[acc,precision,npv,sensitivity,specificity,mcc,roc_auc]
print(result)
row=y_score.shape[0]
y_sparse=utils.to_categorical(y)
yscore_sum = pd.DataFrame(data=y_score)
yscore_sum.to_csv('y_score_RF_Lasso1.csv')
ytest_sum = pd.DataFrame(data=y_sparse)
ytest_sum.to_csv('y_test_RF_Lasso1.csv')
fpr, tpr, _ = roc_curve(y_sparse[:,0], y_score[:,0])
auc_score=result[6]
lw=2
plt.plot(fpr, tpr, color='darkorange',
lw=lw, label='RF ROC (area = %0.4f%%)' % auc_score)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.05])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig(os.path.join('/home/ubuntu/Documents/_Loov1.png'), dpi=600, bbox_inches='tight')
plt.show()
data_csv = pd.DataFrame(data=result)
data_csv.to_csv('RF_Lasso1.csv')