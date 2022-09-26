#!/usr/bin/env python3

import sklearn
import pandas as pd
import numpy as np
import pickle
#import matplotlib.pyplot as plt
from scipy.stats import sem

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.svm import LinearSVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import VotingRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import StackingRegressor
from sklearn.decomposition import PCA

def run_model(model_name, x_scaled, y_train, x_test_scaled, y_test, print_results = True, make_plot = False,prefix ="FF"):
    if print_results == True:
        print(model_name)

    model_name.fit(x_scaled,y_train)
    y_prediction_test = model_name.predict(x_test_scaled)
    y_prediction_train = model_name.predict(x_scaled)
    ereg=model_name
    test_abs_error = []
    train_abs_error = []
    
    if print_results == True:
        print('test ' + str(model_name.score(x_test_scaled,y_test.ravel())))
        print('train '+str(model_name.score(x_scaled,y_train.ravel())))

        for i in range(len(y_test)):
            test_abs_error.append(abs(y_prediction_test[i]-y_test[i]))
        print('test ' + str(np.median(test_abs_error)))

        for i in range(len(y_train)):
            train_abs_error.append(abs(y_prediction_train[i]-y_train[i]))
        print('train '+str(np.median(train_abs_error)))
    if make_plot == True:
        plot_predict_real(y_prediction_train, y_train, train_abs_error, ereg, "{}.Train_statistics".format(prefix))
        plot_predict_real(y_prediction_test, y_test, test_abs_error, ereg, "{}.Test_statistics".format(prefix))
        plot_abs_error(y_prediction_train, y_train, "{}.Absolute_error_test".format(prefix))
        plot_abs_error(y_prediction_train, y_train, "{}.Absolute_error_train".format(prefix))
       
    return model_name

def plot_predict_real(y_predict, y_actual, abs_error_list, model_score, figure_name):
    plt.figure()
    fig1, ax = plt.subplots(figsize = (9, 9))
    ax.scatter(y_predict, y_actual, s=60, alpha=0.7, edgecolors="k")
    ax.annotate("r = {}\nMean abs error: {} +/- {} \n R2: {}".format(np.corrcoef(y_predict, y_actual)[0,1],
                 np.median(abs_error_list), sem(abs_error_list), model_score), (7,15))
    ax.plot(np.unique(y_predict), np.poly1d(np.polyfit(y_predict, y_actual, 1))(np.unique(y_predict)), color='red')
    plt.xlabel("Predicted fetal fraction")
    plt.ylabel("Real FFY value")
    plt.savefig(figure_name+'.png')
    
def plot_abs_error(y_predict, y_actual, figure_name):
    errors = []
    for i in range(len(y_predict)):
        errors.append(y_predict[i]-y_actual[i])
    plt.figure()
    fig2, ax2 = plt.subplots(figsize = (9, 9))
    ax2.hist(errors, bins = 25)
    plt.xlabel("Predicted FFY - Actual FFY")
    plt.savefig(figure_name+'.png')
    


def train_model(training_input_df, model_list, model_names, train_test_proportion, max_FFY, pca_components, print_results = True,prefix="FF"):

    model_ffy=False
    if not "FFY" in training_input_df:
       df_temp = training_input_df.transpose()
       df = df_temp[df_temp["FFY"] < max_FFY]
       df = df[df.FFY != 0.0]
    else:
       model_ffy=True
       df_temp=training_input_df
       df = df_temp[df_temp["FFY"] < max_FFY]

    predictors = list(set(list(df.columns))-set(["FFY","Individual"]))

    #divide data into train and test set
    n_rows = int(len(df.index)*train_test_proportion)
    train = df.iloc[:n_rows,:]
    test = df.iloc[n_rows:,:]
    x_train = train[predictors].values
    x_test = test[predictors].values
    y_train = train['FFY'].values
    y_test = test['FFY'].values
 
    #scale predictor data
    if not model_ffy:
       scaler = StandardScaler().fit(x_train)
       x_scaled = scaler.transform(x_train)
       x_test_scaled = scaler.transform(x_test)
       pca = PCA(n_components=pca_components)
       pca.fit(x_scaled)
       x_scaled = pca.transform(x_scaled)
       x_test_scaled = pca.transform(x_test_scaled)


    else:
       scaler = StandardScaler().fit(df_temp[predictors].values )
       pca = PCA(n_components=pca_components)
       pca.fit(scaler.transform(df_temp[predictors].values) )

       x_scaled = scaler.transform(x_train)
       x_test_scaled = scaler.transform(x_test)
       x_scaled = pca.transform(x_scaled)
       x_test_scaled = pca.transform(x_test_scaled)

    print(len(x_scaled))

    estimators = []
    if len(model_list) > 1:
        for i in range(len(model_list)):
            run_model(model_list[i], x_scaled, y_train, x_test_scaled, y_test,prefix)
            estimators.append(("reg{}".format(i),run_model(model_list[i], x_scaled, y_train, x_test_scaled, y_test, print_results)))

        ereg = VotingRegressor(estimators=estimators)
        if print_results == True:
            print(ereg)
        ereg.fit(x_scaled,y_train)
        y_prediction_test = ereg.predict(x_test_scaled)
        y_prediction_train = ereg.predict(x_scaled)

        test_abs_error = []
        train_abs_error = []
        for i in range(len(y_test)):
            test_abs_error.append(abs(y_prediction_test[i]-y_test[i]))
        for i in range(len(y_train)):
            train_abs_error.append(abs(y_prediction_train[i]-y_train[i]))

        test_score = ereg.score(x_test_scaled,y_test.ravel())
        train_score = ereg.score(x_scaled,y_train.ravel())

        if print_results == True:
            print('test ' + str(ereg.score(x_test_scaled,y_test.ravel())))
            print('train '+str(ereg.score(x_scaled,y_train.ravel())))
            print('test ' + str(np.median(test_abs_error)))
            print('train '+str(np.median(train_abs_error)))
        plot_predict_real(y_prediction_train, y_train, train_abs_error, train_score, "{}_Train_statistics".format(prefix))
        plot_predict_real(y_prediction_test, y_test, test_abs_error, test_score, "{}_Test_statistics".format(prefix))
        plot_abs_error(y_prediction_train, y_train, "{}_Absolute_error_test".format(prefix))
        plot_abs_error(y_prediction_train, y_train, "{}_Absolute_error_train".format(prefix))

        with open("{}_model.pkl".format(prefix), 'wb') as file:
            pickle.dump({'model':ereg, 'pca':pca, 'scaler':scaler, 'predictors':predictors}, file)
    else:
        model = run_model(model_list[0],x_scaled, y_train, x_test_scaled, y_test, print_results,prefix)
        print(model)
        with open("{}_model.pkl".format(prefix), 'wb') as file:
            pickle.dump({'model':model, 'pca':pca, 'scaler':scaler, 'predictors':predictors}, file)

    
  

    

#adam = MLPRegressor(hidden_layer_sizes=(5,5),solver='sgd',max_iter=10000,activation="logistic",random_state=40,n_iter_no_change=5)
#LR = LinearRegression()
#lvsm = LinearSVR(random_state=0, tol=1e-5,C=10,max_iter=1000000)
#randomforest=RandomForestRegressor(max_depth=10, random_state=0,n_estimators=100)

#train_model('/proj/sens2017106/esmee2/combined_bins.csv', [adam,LR,lvsm,randomforest], ['adam','LR','lvsm','randomforest'], 0.8, 20, 10, print_results = True)
