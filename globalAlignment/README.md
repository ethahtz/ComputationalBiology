Readme file for align.py
Author: Etha Hua
Date: Feb 8th, 2022
Purpose: Implementation of the global alignment algorithm
with a naive scoring function and a linear gap genalty.

Usage: 

For default scoring scheme (M=4, m=-2, g=-2):
python3 align.py < [inputFileName]
e.g.: python3 align.py < aligntest.input1

For useridefined scoring scheme:
python3 align.py [M] [n] [g] < [inputFileName]
e.g.: python3 align.py 6 2 -2 < aligntest.input1
