from MotifFinder import *
import math

path = r"Data/MedianString/inputs/dataset_30304_9.txt"

f = open(path)
lines = f.readlines()
k = int(lines[0].rstrip())

dnas = []
for line in lines[1:]:
    dnas.append(line.rstrip())

print(dnas)
print(Median_String(dnas, k))


f.close()

# dnas = ["AAATTGACGCAT", "GACGACCACGTT",
#          "CGTCAGCGCCTG", "GCTGAGCACCGG",
#          "AGTTCGGGACAG"]

# print(Median_String(dnas, 3))

# pattern = "AAA"
#
# # dnas = ["AAAAA", "AAAAA", "AAAAA"]
# k= 5
# d = 2

# path = r"Data/MotifEnumeration/outputs/output_7.txt"
# out_path = "output.txt"

# f = open(path)

# data = f.read().rstrip().split(' ')

# print(len(data))

# f.close()

# with open(out_path, 'w') as o:
#     results = sorted(Motif_Enumeration(dnas, k, d))
#     for res in results:
#         o.write(f"{res} ")


