from nolearn.dbn import DBN
import numpy as np
from sklearn.datasets import fetch_mldata
from gdbn.dbn import buildDBN
from gdbn import activationFunctions
from sklearn.preprocessing import MinMaxScaler
from sklearn.cross_validation import KFold
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix, precision_score, recall_score
import csv
import sys
from time import gmtime, strftime


with open('GenotipeBiner.csv', 'r') as g:
    read = csv.reader(g)
    geno = np.array(list(read), dtype='float64')
    g.close()

with open('Phenotype RUB Binary 2 Class.csv', 'rb') as f:
    read = csv.reader(f)
    feno = np.array(list(read))
    f.close()

basic_x = geno/3.0
basic_y = feno.flatten()
basic_y.astype("int0")

dbn_model = DBN([basic_x.shape[1], 100, 100, 100, 2],
                    dropouts=0.01,
                    output_act_funct=activationFunctions.Softmax(),
                    learn_rates=0.01,
                    learn_rates_pretrain=0.001,
                    minibatch_size=9,
                    learn_rate_decays=0.9,
                    learn_rate_minimums=0.0001,
                    epochs_pretrain=5000,
                    epochs=5000,
#                    momentum= self.momentum,
                    real_valued_vis=True,
                    use_re_lu=True,
                    verbose=1)



akurasi = 0.0
#presisi = 0.0
#recall = 0.0
n = 0
waktu_eksekusi = strftime("%Y-%m-%d %H.%M.%S", gmtime())
filename = 'output-' + waktu_eksekusi + 'shuffleFalse.txt'

for i in range(2):
    kf = KFold(81, n_folds=9, shuffle=False)
    a = 0
    for train_index, test_index in kf:
        a += 1
        X_train, X_test = basic_x[train_index], basic_x[test_index]
        y_train, y_test = basic_y[train_index], basic_y[test_index]
        print 'Building model iterasi ke- %d fold ke- %d.\n' % (i, a)
        dbn_model.fit(X_train, y_train)

        y_true, y_pred = y_test, dbn_model.predict(X_test) # Get our predictions
        orig_stdout = sys.stdout
        sys.stdout = open(filename, 'a+')
        print '============================================================================================='
        #print(classification_report(y_true, y_pred)) # Classification on each digit

        print 'Output sebenarnya:', y_true
        print '\nOutput prediksi:', y_pred

        print 'Confusion matrix:\n',(confusion_matrix(y_true, y_pred))
        print '\n'
        print(classification_report(y_true, y_pred))
        print '\nAkurasi:', accuracy_score(y_true, y_pred)
        #print '\nPresisi:', precision_score(y_true, y_pred)
        #print '\nRecall:', recall_score(y_true, y_pred)
        print '============================================================================================='
        akurasi += accuracy_score(y_true, y_pred)
        #presisi += precision_score(y_true, y_pred)
        #recall += recall_score(y_true, y_pred)
        n += 1
        sys.stdout.close()
        sys.stdout = orig_stdout
        print '\nAkurasi:', accuracy_score(y_true, y_pred)
        print strftime("%Y-%m-%d %H:%M:%S", gmtime())

orig_stdout = sys.stdout
sys.stdout = open(filename, 'a+')
print 'Jumlah percobaan sebanyak:', n
print '\nRata-rata akurasi adalah:', akurasi/n
waktu_selesai = strftime("%Y-%m-%d %H.%M.%S", gmtime())
print '\nWaktu start:', waktu_eksekusi
print '\nWaktu finish:', waktu_selesai
#print '\nRata-rata presisi adalah:', presisi/n
#print '\nRata-rata recall adalah:', recall/n
sys.stdout.close()
sys.stdout = orig_stdout


print 'SELESAI\nJumlah percobaan sebanyak:', n
print '\nRata-rata akurasi adalah:', akurasi/n