Readme file for deBruijnGraph.py  
Author: Etha Hua  
Date: March 14th, 2022  
Purpose: Implementation of the de Bruijn graph for sequence assembly  
  
Usage:  
For producing good reads (mode "g"):  
  
python3 deBruijnGraph.py [sequence reads filename] g [k-value]  
  
e.g. python3 deBruijnGraph.py sequence_reads g 31  
  
For finding contigs (mode"c"):  
  
python3 deBruijnGraph.py [good reads filename] c [k-value] (min length contig)  
  
e.g. python3 deBruijnGraph.py good_reads c 31 100  
