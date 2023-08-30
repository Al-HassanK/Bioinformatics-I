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
            
                    

