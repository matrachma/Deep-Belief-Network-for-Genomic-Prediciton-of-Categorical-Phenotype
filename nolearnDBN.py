from nolearn.dbn import DBN
import numpy as np
from sklearn.datasets import fetch_mldata
from gdbn.dbn import buildDBN
from gdbn import activationFunctions
from sklearn.preprocessing import MinMaxScaler
from sklearn.cross_validation import train_test_split
from sklearn.metrics import classification_report, accuracy_score
import csv


def train_test_prep():

    with open('GenotipeBiner.csv', 'r') as g:
        read = csv.reader(g)
        geno = np.array(list(read), dtype='float64')
        g.close()

    with open('Phenotype RUB Binary 2 Class.csv', 'rb') as f:
        read = csv.reader(f)
        feno = np.array(list(read))
        f.close()

    basic_x = geno
    basic_y = feno.flatten()

    #min_max_scaler = MinMaxScaler() # Create the MinMax object.
    #basic_x = min_max_scaler.fit_transform(basic_x) # Scale pixel intensities only.


    #x_train, x_test, y_train, y_test = basic_x, basic_x, basic_y, basic_y
    x_train, x_test, y_train, y_test = train_test_split(basic_x/2.0, basic_y.astype("int0"),
                            test_size = 0.1, random_state = 0) # Split training/test.

    #return x_train, x_test, y_train, y_test
    return x_train, x_test, y_train, y_test, basic_x, basic_y.astype("int0")

x_train, x_test, y_train, y_test, x, y = train_test_prep()

#dbn_model = DBN([x_train.shape[1], 1500, 1500, 2],
#                learn_rates = 0.01,
#                learn_rate_decays = 0.9,
#                epochs = 1000,
#                verbose = 3)

dbn_model = DBN([x_train.shape[1], 5000, 2500, 1250, 500],
                    #dropouts=0.01,
                    output_act_funct=activationFunctions.Sigmoid(),
                    learn_rates=0.01,
                    learn_rates_pretrain=0.001,
#                    minibatch_size=9,
#                    learn_rate_decays=0.9,
#                    learn_rate_minimums=0.0001,
                    epochs_pretrain=500,
                    epochs=500,
#                    momentum= self.momentum,
#                    real_valued_vis=True,
#                    use_re_lu=True,
                    verbose=2)


dbn_model.fit(x_train, y_train)

y_true, y_pred = y_test, dbn_model.predict(x_test) # Get our predictions
print(classification_report(y_true, y_pred)) # Classification on each digit
print y_true
print y_pred

print '\nThe accuracy is:', accuracy_score(y_true, y_pred)
