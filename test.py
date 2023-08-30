from MotifFinder import *

dnas = ["CGTTACAGATGCTAGGAAAGTAGCT", "AGCTCCACGAGAATGTCTGTCAGAA",
         "CGAATTGGATCGCTCAAACTGAGAG", "TCGTGATTTCCATCTGAGGGAAGCT",
         "TATTGGGATTTAAGTGAAAGCCGCT", "AAGTACAATGTTCCCGAGGGGTACG"]

# dnas = ["AAAAA", "AAAAA", "AAAAA"]
k= 5
d = 2

path = r"Data/MotifEnumeration/outputs/output_7.txt"
out_path = "output.txt"

f = open(path)

data = f.read().rstrip().split(' ')

print(len(data))

f.close()

with open(out_path, 'w') as o:
    results = sorted(Motif_Enumeration(dnas, k, d))
    for res in results:
        o.write(f"{res} ")

