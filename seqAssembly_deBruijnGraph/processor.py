import numpy as np

fstream = open("contig_lengths", "r")

contigLen = []

for line in fstream:
    contigLen.append(int(line))

contigLen.sort(reverse=True)

# print(contigLen)


s2 = 0

n50 = np.sum(contigLen)/2

for i in range(len(contigLen)):
    s2 += contigLen[i]
    if s2 >= n50:
        print("the ", i, "exceeds", n50, "total sum:", n50 * 2)
        print("sum of contig lengths up to this:",s2,",N-50 size:", contigLen[i])
        break

print("Max contig Length:", contigLen[0])

# numShort = 0

# for elem in contigLen:
#     if elem < 100:
#         numShort += 1
# print(numShort/len(contigLen))