from Methods import *
import collections
import numpy
import random

### Pre-conditions: Takes a dna string, integer k_mer which is the length of the motif...
### Post-conditions: Returns a list of all k_mer patterns in the dna strign
def K_mer_Divider(dna_string, k_mer):
    k_mers = []
    for i in range(0, len(dna_string) - k_mer + 1):
        k_mers.append(dna_string[i:i+k_mer])
    
    return k_mers

### Pre-conditions: Takes a dna pattern, alist of dna sequences, d (number of mismatches)...
### Post-conditions: Returns true if the pattern exists in all dna strings with at most d mismtaches... 
def Check_Pattern_Existence(pattern, dna_strings, d):
    flags = []
    for dna_string in dna_strings:
        if Approximate_Matched_Positions(pattern, dna_string, d):
            flags.append(True)
        else:
            flags.append(False)
    
    if False in flags:
        return False
    else:
        return True

### Pre-conditions: Takes a collection (list) of dna strings, two integers k_mer which is the length of motif, and d the maximum
# number of mismatches...
### Post-conditions: Returns all (k_mer, d)-motifs in the dna collection... 
def Motif_Enumeration(dna_strings, k_mer, d):
    """
        This function is Exhaustive alogrithm to find a motif of length k_mer if it appears in every DNA string, of a passed
        collection of strings, with at most d mismatches.

    """
    patterns = []
    # Generate all k-mer patterns in the first dna string
    first_string_k_mers = K_mer_Divider(dna_strings[0], k_mer)
    
    for pattern in first_string_k_mers:
        for neighbor in Neighbors(pattern, d): # Generate all neighbors for each k-mer pattern
            if len(neighbor) == k_mer:
                # Check whether each neighbor for that pattern appears in each string from dna strings with at most d mismatches
                if Check_Pattern_Existence(neighbor, dna_strings, d):
                    patterns.append(neighbor)

    patterns = set(patterns) # remove duplicates

    return list(patterns) # returns a list of all disinct (k-mer, d)-motifs
            
### Pre-conditions: a dna pattern and a list of dna strings...
### Post-conditions: The sum of minimum distances between pattern and all dna strings, in other words the minimum distance
# between a given pattern and collection of DNA strings...                   
def Distance_Between_Pattern_And_Strings(pattern, dna_strings):
    k = len(pattern) 
    distance = 0
    for dna_string in dna_strings:
        hamming_distance = float("inf")
        # The following loop do this: d(pattern, Text) = min(all k-mers pattern` in Text) Hamming_Distance(patter, pattern`)
        # In other words it finds the minimum hamming distance value between the given pattern and a substring in Text
        for k_mer in K_mer_Divider(dna_string, k):
            if hamming_distance > Hamming_Distance(pattern, k_mer):
                hamming_distance = Hamming_Distance(pattern, k_mer)
        
        distance += hamming_distance

    return distance

### Pre-conditions: Takes an integer k that represent the length of the DNA pattern...
### Post-conditions: Returns a list of all possible DNA sequences of length k...
def All_Strings(k):
    start_pattern = k * 'A'

    return Neighbors(start_pattern, k)

### Pre-conditions: Takes a list of DNA strings, and an integer k...
### Post-conditions: Returns a k-mer pattern that minimizes d(pattern, list_of_dna_strings) among all possible choices of k-mers...
def Median_String(dna_strings, k):
    distance = float("inf")
    median = ""

    patterns = All_Strings(k) # Find all possible patterns of size k.
    for pattern in patterns:
        if distance > Distance_Between_Pattern_And_Strings(pattern, dna_strings):
            distance = Distance_Between_Pattern_And_Strings(pattern, dna_strings)
            median = pattern

    return median

### Pre-conditions: A DNA string, an integer k, and 4 x k matrix profile...
### Post-conditions: A profile-most probable k-mer in DNA string...
def Find_Profile_Most_Probable_K_mer(dna_string, k, profile_matrix):
    all_k_mers = K_mer_Divider(dna_string, k)

    profile_matrix_indexder = {'A':0, 'C':1, 'G':2, 'T':3}
    k_mers_probabilities = []

    for k_mer in all_k_mers:
        probability = 1.0
        for i in range(k):
            row_index = profile_matrix_indexder[k_mer[i]]
            probability *= profile_matrix[row_index][i]
        
        k_mers_probabilities.append(probability)
    
    most_probable_k_mer_index = k_mers_probabilities.index(max(k_mers_probabilities))

    return all_k_mers[most_probable_k_mer_index]

### Pre-conditions: takes a list of DNA motifs...
### Post-conditions: Returns a numpy array of size 4 * k (the length of one motif) which contains the occurrences of each
# nucleotide in each column of the given motifs list...
def Count(motifs):
    """
        In this function Laplace's Rule of Succesion is applied 
    """

    k = len(motifs[0]) # The length of motif
    t = len(motifs) # the number of motifs in the list...

    count_matrix = numpy.zeros((4, k)) # will hold the final result...
    count_matrix_indexder = {'A':0, 'C':1, 'G':2, 'T':3}

    for i in range(k):
        column_sequence = ""
        for j in range(t):
            column_sequence += motifs[j][i] # holds each column sequence of the motifs list
        
        nucleotide_counter = dict(collections.Counter(column_sequence)) # calculating nucleotide occurences in the column...
        # assign the values to the count_matrix...
        for key in count_matrix_indexder.keys():
            if not nucleotide_counter.get(key, None):
                count_matrix[count_matrix_indexder[key], i] = 1
            else:
                 count_matrix[count_matrix_indexder[key], i] = nucleotide_counter[key] + 1

    return count_matrix
### Pre-conditions: takes a list of DNA motifs...
### Post-conditions: Returns a list of size 4 * k (the length of one motif) which contains the occurrences of each
# nucleotide in each column of the given motifs list divided by the number of motifs...
def Profile(motifs):
    k = len(motifs[0]) # The length of motif
    t = len(motifs) # the number of motifs in the list...

    count_matrix = Count(motifs) # Generate Count(motifs) matrix
    profile_matrix = [] # will hold the final results...

    for i in range(count_matrix.shape[0]):
        row_profile = [] # the row probabilities...
        for j in range(count_matrix.shape[1]):
            row_profile.append(count_matrix[i, j] / float(t))
        
        profile_matrix.append(row_profile)

    return profile_matrix

### Pre-conditions: takes a list of DNA motifs...
### Post-conditions: Retruns the sum of unpopular letters in the motifs matrix...
def Score(motifs):
    count_matrix = Count(motifs)
    t = len(motifs)
    return numpy.sum([(t - max_value) for max_value in numpy.max(count_matrix, axis=0)])
    
### Pre-conditions: Takes a list of DNA strings, and integers k and t...
### Post-conditions: Returns a collection of strings best_motifs...
def Greedy_Motif_Search(dna_strings, k, t):
    best_motifs = [K_mer_Divider(dna_string, k)[0] for dna_string in dna_strings]

    for k_mer in K_mer_Divider(dna_strings[0], k):
        motifs = [k_mer]
        for i in range(1, t):
            profile = Profile(motifs)
            motifs.append(Find_Profile_Most_Probable_K_mer(dna_strings[i], k, profile))

        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs

    return best_motifs

### Pre-conditions: Takes a DNA sequence and integer k...
### Post-conditions: Returns a random sequence of length k from that sequence...
def Generate_Random_K_mers(dna_strings, k):
    randMotifs = []
    t = len(dna_strings)

    for i in range(t):
        x = random.randint(0, t)
        randMotifs.append(dna_strings[i][x:x+k])

    return randMotifs

### Pre-conditions: Takes a profile matrix, dna strings and an integer k...
### Post-conditions: Returns list of k-mers formed by the Profile-most probable k-mers in each string from dna strings...
def Motifs(profile_matrix, dna_strings, k):
    motifs = []
    for dna_string in dna_strings:
        motifs.append(Find_Profile_Most_Probable_K_mer(dna_string, k, profile_matrix))

    return motifs

### Pre-conditions: Integers k and t, and a list of strings Dna...
### Post-conditions: A list of bestMotifs...
def Randomized_Motif_Search(dna_strings, k):
    motifs = Generate_Random_K_mers(dna_strings, k)
    best_motifs = motifs

    while True:
        profile = Profile(motifs)
        motifs  = Motifs(profile, dna_strings, k)

        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs
        

def Run_Randomized_Motif_Search(dna_strings, k, t, num_of_iterations=1000):
    motifs = Randomized_Motif_Search(dna_strings, k)
    best_motifs = motifs

    for i in range(num_of_iterations):
        motifs = Randomized_Motif_Search(dna_strings, k)
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        else:
            final_best_motifs = best_motifs
        
    return final_best_motifs