Readme file for gibbsSampler.py 
Author: Etha Hua
Date: Feb 24th, 2022
Purpose: Implementation of the Gibbs Sampler algorithm


Usage for gibbsSampler-i.py:

For default motif length(6):
`python3 gibbsSampler-i.py < [inputFileName]`
e.g.: `python3 gibbsSampler-i.py < Gibbs.fasta`

For user-defined motif length:
`python3 gibbsSampler-i.py [motifLength] < [inputFileName]`
e.g.: `python3 gibbsSampler-i.py 10 < Gibbs.fasta`

For user-defined motif length and convergence criterion:
`python3 gibbsSampler-i.py [motifLength] [conv] < [inputFileName]`
e.g.: `python3 gibbsSampler-i.py 10 50 < Gibbs.fasta`
