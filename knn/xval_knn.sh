#!/bin/sh
python3 knn-classifier.py GSE994-train.txt GSE994-test.txt 1 xv > Prob5-1NNxval.txt
python3 knn-classifier.py GSE994-train.txt GSE994-test.txt 3 xv > Prob5-3NNxval.txt