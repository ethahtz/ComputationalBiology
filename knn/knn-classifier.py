"""
knn-classifier.py

Description: Implementation of the k-nearest-neighbours
classifier for smoking behavior given features of 
gene expression levels. Supports prediction of the 
testing samples as well as cross-validation within
the training samples
    
Created by Etha Hua, April 28, 2022
"""
from cgi import test
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy import stats
import sys

# get_data
# Purpose: Extract training and testing samples from files
#          and return them in numpy arrays
# Parameters: train_fname - file name of the file with training samples
#             test_fname - file name of the file with testing samples
# Returns: train_features - numpy array of the features of the training
#                           samples
#          train_labels - numpy array of the labels of the training
#                         samples
#          test_features - numpy array of the features of the testing 
#                          samples 
def get_data(train_fname, test_fname):
    X = pd.read_csv(train_fname, delimiter='\t')
    XT = pd.read_csv(test_fname, delimiter='\t')

    train_data = X.values.T
    train_features = np.array(train_data[:,:-1], dtype=float)
    train_labels = train_data[:, -1]

    test_data = XT.values.T
    test_features = np.array(test_data[:,:-1], dtype=float)

    return train_features, train_labels, test_features

# knn_decider
# Purpose: Given an array of distances (from a testing sample to all 
#          training samples) and an array of labels of the training
#          samples, decide whether that testing samples belongs to
#          which label according to the k_value and majority vote 
# Parameters: dists - an array of float numbers representing the 
#                     distances from a testing sample to all 
#                     training samples 
#             labels - an array of labels of the training samples
#             k_value - the number of neighbours one set for the 
#                       classifier
# Returns: the decided label for the testing sample
def knn_decider(dists, labels, k_value):
    knns = []
    for i in range(k_value):
        knns.append((dists[i], labels[i]))
        knns.sort()

    for i in range(k_value,len(dists)):
        index = k_value
        while index >= 1 and dists[i] < knns[index-1][0]:
            index -= 1
        if index != k_value:
            knns.insert(index, (dists[i], labels[i]))
            knns.pop(k_value)

    labels = []
    for pairs in knns:
        labels.append(pairs[1])
        
    return stats.mode(labels)[0][0]

# knn_predict
# Purpose: Predict an array of testing samples' labels given
#          some labeled training samples and a k-value 
# Parameters: train_features - features of the training samples
#             train_labels - labels of the training samples
#             test_features - features of the testing samples
#             k-value - hyper parameter of num of neighbours to 
#                       check
# Returns: an array of labels for the testing samples 
def knn_predict(train_features, train_labels, test_features, k_value):
    all_dists = distance.cdist(test_features, train_features, 'euclidean')
    predictions = []
    for i in range(len(test_features)):
        predictions.append(knn_decider(all_dists[i],train_labels,k_value))
    return predictions

# x_validation
# Purpose: Cross-validate the algorithm given a set of labeled samples
# Parameters: features - features of the labeled samples 
#             labels - labels of the labeled samples
#             k-value - hyper parameter of num of neighbours to 
#                       check
# Returns: an array of labels for the samples
def x_validation(features, labels, k_value):
    fold_size = 6
    predicted = []
    numItr = int(len(features) / fold_size)
    for i in range(numItr):
        currTrain = np.delete(features, range(i * fold_size, i * fold_size + fold_size), axis=0)
        currLabels = np.delete(labels, range(i * fold_size, i * fold_size + fold_size))
        currTest = np.take(features, range(i * fold_size, i * fold_size + fold_size), axis=0)
        curr_prediction = knn_predict(currTrain, currLabels, currTest, k_value)
        predicted.append(curr_prediction)
    result = np.array(predicted).reshape(1,-1)[0]
    return result

# calc_accuracy
# Purpose: Calculate the accuracy of prediction given ground truth
# Parameters: predicted - an array of predicted labels
#             true-labels - an array of true labels
# Returns: a float number between 0 and 1 representing the 
#          accuracy of the prediction
def calc_accuracy(predicted, true_labels):
    numCorrect = 0
    for i in range(len(predicted)):
        if predicted[i] == true_labels[i]:
            numCorrect += 1
    return float(numCorrect/len(predicted))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 knn-classifier.py [TRAIN.txt] [TEST.txt] [K-VALUE] [test/xv]")
        exit(1)
    train_fname, test_fname, k_val = sys.argv[1], sys.argv[2], int(sys.argv[3])
    train_features, train_labels, test_features = get_data(train_fname, test_fname)
    modeProgram = sys.argv[4]
    if modeProgram == "test":
        predicted_labels = knn_predict(train_features, train_labels, test_features, k_val)
        test_headers = pd.read_csv(sys.argv[2], nrows=1, delimiter='\t').columns.to_list()
        for i in range(len(test_features)):
            print(test_headers[i], predicted_labels[i])
    elif modeProgram == "xv":
        vldrst = x_validation(train_features, train_labels, k_val)
        test_headers_vld = pd.read_csv(sys.argv[1], nrows=1, delimiter='\t').columns.to_list()
        for i in range(len(train_labels)):
            print(test_headers_vld[i], vldrst[i])
        # Uncomment the following line for printing the accuracy rate of
        # Cross-validation predicted labels corresponding to the features
        # of the training samples 

        # print(calc_accuracy(vldrst, train_labels))
    else:
        print("Unsupported Mode, please use 'test' for prediction or 'xv' for cross-validation")