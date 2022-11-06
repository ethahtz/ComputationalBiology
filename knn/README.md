Readme file for HW5 Part 2  
Author: Etha Hua  
Date: April 28th, 2022  
Purpose: Implementation of a K nearest neighbours classifier 
for smoking behavior identification based on gene expressions

-----Running Instructions-----  
Run the KNN classifier with:  
  
python3 knn-classifier.py [TRAIN.txt] [TEST.txt] [K-VALUE] [test/xv]  
  
where   
  
TRAIN.txt is replaced by the filename of the training samples,  
TEST.txt is replaced by the filename of the testing samples,  
K-VALUE is replaced by the k-values - number of neighbours for the classifier,  
and [test] for predicting the labels of the testing samples or 
[xv] for cross-validation within the training samples.

For producing requested outputfiles, the user could run the 
shell scripts written for that:
  
For producing 1-NN and 3-NN predictions on the testing samples:

sh predict_test.sh
  
For producing 1-NN and 3-NN cross-validation predictions:

sh xval_knn.sh
