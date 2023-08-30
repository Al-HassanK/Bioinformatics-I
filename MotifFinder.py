from Methods import *

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
            
### Pre-conditions: a dna patterns and a list of dna strings...
### Post-conditions: The sum of minimum distances between pattern and all dna strings...                   
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


    