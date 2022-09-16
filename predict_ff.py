
import pandas as pd
import sklearn
import pickle

def predict_FF(df):
    ind_list = df.columns[1:]
    with open("FF_prediction_model.pkl", 'rb') as file:
        data = pickle.load(file)
    pca = data["pca"]

    model=data["model"]
    scaler = data["scaler"]
    predictors = data["predictors"]      

    print("Individual,FFY,FFpredicted)   
    for ind in ind_list:
        test = df.[ind]
        FFY = df.[ind]['FFY']
        x_test = test[predictors].values
        x_scaled = scaler.transform(x_test.reshape(1,-1))
        x_scaled = pca.transform(x_scaled)
        FF_predict = model.predict(x_scaled)
        print(ind + ', '  + str(FFY)+', ' + str(FF_predict[0])+'\n')    

