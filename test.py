from MotifFinder import *
import math
import scipy.stats

# motifs = ["TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC"
#           ,"TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"]

# motifs = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", 
#           "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
# k = 8
# t = 5

# res = Run_Randomized_Motif_Search(motifs, k, t)
# for r in res:
#     print(r, end=' ')

# dna = []
# best = Greedy_Motif_Search(motifs, k, t)
# print(best)


path = r"Data/RandomizedMotifSearch/inputs/dataset_30307_5.txt"
f = open(path)
lines = f.readlines()
k, t = lines[0].rstrip().split(' ')
profile = []
dna_strings = []
for line in lines[1:]:
    for dna_string in line.rstrip().split(' '):
        dna_strings.append(dna_string)

res = Run_Randomized_Motif_Search(dna_strings, int(k), int(t))
for r in res:
    print(r, end=' ')    
# print(res)

# print(dna_strings[0][25:])

# for res in Greedy_Motif_Search(dna_strings, int(k), int(t)):
#     print(res, end=' ')
# print(Find_Profile_Most_Probable_K_mer(dna, k, profile))
# # print(profile)
# f.close()

# p = [
#     [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
#     [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
#     [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
#     [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
# ]
# print(scipy.stats.entropy(p[0]))

# path = r"Data/MedianString/inputs/dataset_30304_9.txt"

# f = open(path)
# lines = f.readlines()
# k = int(lines[0].rstrip())

# dnas = []
# for line in lines[1:]:
#     dnas.append(line.rstrip())

# print(dnas)
# print(Median_String(dnas, k))


# f.close()

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


# def Entropy(l):
#     entropy = [(p * math.log2(p)) for p in l]
#     return (- sum(entropy)) 

# print(Entropy([0.5, 0.0, 0.0, 0.5]))
# print(Entropy([0.25, 0.25, 0.25, 0.25]))
# print(Entropy([0.0, 0.0, 0.0, 1]))
# print(Entropy([0.25, 0.0, 0.5, 0.25]))