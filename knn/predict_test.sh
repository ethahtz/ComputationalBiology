#!/bin/sh
python3 knn-classifier.py GSE994-train.txt GSE994-test.txt 1 test > Prob5-1NNoutput.txt
python3 knn-classifier.py GSE994-train.txt GSE994-test.txt 3 test > Prob5-3NNoutput.txt