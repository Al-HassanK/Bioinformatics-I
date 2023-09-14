"""
    This File contains the functions in the interactive text for weeks 1 and 2. The methods and algorithms defined here are used to solve popular string problems
    like, the finding clumps problem and approximate matching. You can use these method to develop a program that takes a bacterial genome then it can finds and locates
    the origin of replication "Ori"...

"""

### Pre-conditions: Takes a text and a pattern you want to look for in the text...
### Post-conditions: Returns the number of occurences of the pattern in the text...
def Pattern_Count(ref_text, pattern):
    count = 0
    
    for i in range(0, (len(ref_text)-len(pattern)) + 1):
        if ref_text[i:len(pattern) + i] == pattern:
            count += 1
    
    return count

### Pre-conditions: Takes a text and a length of k_mer...
### Post-conditions: Returns a table "Dictionary" of each pattern of size k_mer in the text with its corresponding count... 
def Frequency_Table(ref_text, k_mer):
    frequency_map = {}
    n = len(ref_text)

    for i in range(0, (n-k_mer) + 1):
        pattern = ref_text[i:k_mer+i] # Get a substring from the ref_text of size k_mer...
        if frequency_map.get(pattern) == None: # Check whether this pattern in the frequency table or not...
            frequency_map[pattern] = 1
        else:
            frequency_map[pattern] += 1
    
    return frequency_map

### Pre-conditions: Takes frequency table returned from Frequency_Table function...
### Post-conditions: Returns the most frequent pattern in the frequency table... 
def Max_Map(frequency_table):
    return max(frequency_table.values())

### Pre-conditions: Takes a text and a length of k_mer...
### Post-conditions: Returns a list of the most frequent k_mer patterns in the ref_text... 
def Frequent_Words(ref_text, k_mer):
    frequent_patterns = []
    freq_map = Frequency_Table(ref_text, k_mer) # Get the frequency table...
    max_value = Max_Map(freq_map) # Get the maximum occurence in the table

    for pattern in freq_map.keys():
        if freq_map[pattern] == max_value:
            frequent_patterns.append(pattern)

    return frequent_patterns

### Pre-conditions: Takes a genome, window_size "Which is the length of the region we will search within", a k_mer, 
# and a threshold t...
### Post-conditions: Returns all distinct patterns of size k_mer forming (window_size, t)-clumps in the genome...
def Find_Clumps(genome, window_size, k_mer, t):
    patterns = set() 
    n = len(genome)

    for i in range(0, (n - window_size) + 1):
        window = genome[i:window_size+i]
        freq_map = Frequency_Table(window, k_mer) # Get the frequency table of that region of the genome
        # Add all the patterns that occurs more than or equals to t to the patterns set... 
        for pattern in freq_map.keys():
            if freq_map[pattern] >= t:
                patterns.add(pattern)
    
    return patterns

### Pre-conditions: Takes a genome and an index i...
### Post-conditions: Returns the skew diagram at the position i...
def Skew_Diagram(genome, i=None):
    """
        This function calculates the difference between the total number of occurrences of G and the total
        number of occurrences of C in the first i nucleotides of the genome, if i was not specified the difference will
        be calculated for the whole genome and if i is set to zero then the skew is just zero...

    """
    if i == None:
        i = len(genome)

    skew_values = [0] # the first value in the skew diagram is always zero...

    # traverse the genome until the i postion while adding the occurrences of Gs and subtracting the occurrences of Cs...
    for index in range(1, i+1):
        if genome[index-1] == 'G':
            skew_values.append(skew_values[index-1] + 1)
        elif genome[index-1] == 'C':
            skew_values.append(skew_values[index-1] - 1)
        else:
            skew_values.append(skew_values[index-1])

    return skew_values

### Pre-conditions: Takes a genome...
### Post-conditions: Retruns all positions i minimizing Skew_Diagrame of that genome among all values of i... 
def Minimun_Skew(genome):
    skew_values = Skew_Diagram(genome) # Build the skew diagram of the genome...
    min_skew = min(skew_values) # Get the minimum value within it...
    positions = []
    # Make a list of all the positions of the most minimum values in the skew diagram...
    for i in range(0, len(skew_values)):
        if skew_values[i] == min_skew:
            positions.append(i)
            break
    
    return positions

### Pre-conditions: Two strings of equal length...
### Post-conditions: The number of mismatches between the two strings...
def Hamming_Distance(seq1, seq2):
    num_of_mismatches = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            num_of_mismatches += 1
    
    return num_of_mismatches

### Pre-conditions: Takes a pattern, a referenece text, and a threshold d...
### Post-conditions: Retruns all starting positions where pattern appears as a substring of ref_text with at most d mismatches...
def Approximate_Matched_Positions(pattern, ref_text, d):
    positions = []
    for i in range(0, len(ref_text)-len(pattern)+1):
        if Hamming_Distance(pattern, ref_text[i:len(pattern)+i]) <= d:
            positions.append(i)

    return positions


### Pre-conditions: Takes a pattern, a referenece text, and a threshold d...
### Post-conditions: Retruns the number of starting positions where pattern appears as a substring of ref_text 
# with at most d mismatches...
def Approximate_Pattern_Count(ref_text, pattern, d):
    return len(Approximate_Matched_Positions(pattern, ref_text, d))


### Pre-conditions: Takes a pttern and a threshold d...
### Post-conditions: Returns the set of all k-mers whose Hamming_Distance from pattern does not exceed d... 
def Neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    
    nucleotides = {'A', 'C', 'G', 'T'}
    neighborhood = set()
    suffix_neighbors = Neighbors(pattern[1:], d)
    for suffix_neighbor in suffix_neighbors:
        if Hamming_Distance(pattern[1:], suffix_neighbor) < d:
            for nucleotide in nucleotides:
                neighborhood.add(nucleotide + suffix_neighbor)
        else:
            neighborhood.add(pattern[0] + suffix_neighbor)
    
    return neighborhood

### Pre-conditions: Takes a reference text, and integers k_mer and d...
### Post-conditions: Returns all most frequent k_mer patterns with up to d mismatches in ref_text...
def Frequent_Words_With_Mismathces(ref_text, k_mer, d):
    patterns = [] # will contain the most frequent k_mer patterns
    freq_map = {} # the frequency table of all k_mer patterns
    n = len(ref_text)

    # This for loop generates all k_mer patterns and their neighbors and populate the frequency table with them...
    for i in range(0, n-k_mer+1):
        pattern = ref_text[i:k_mer+i] 
        neighborhood = Neighbors(pattern, d) # Generate all neighbors for that pattern as well it will contain the original pattern
        for neighbor in neighborhood:
            if freq_map.get(neighbor) == None:
                freq_map[neighbor] = 1
            else:
                freq_map[neighbor] += 1

    m = Max_Map(freq_map) # Get the maximum frequency between all patterns
    # Find all patterns that have a maximum occurence in the text...
    for pattern in freq_map.keys():
        if freq_map[pattern] == m:
            patterns.append(pattern)
    
    return patterns


### Pre-conditions: Takes a DNA sequence...
### Post-conditions: Returns its reverse complement...
def Reverse_Complement(sequence):
    complemets = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    reveresed_sequence = ''
    for n in sequence:
        reveresed_sequence = complemets[n] + reveresed_sequence
    
    return reveresed_sequence


def Max_Map_With_Reverse_Complements(frequency_map):
    combined_values = []
    for key in frequency_map.keys():
        key_reverse_complement = Reverse_Complement(key)
        combined_values.append(frequency_map[key] + frequency_map.get(key_reverse_complement, 0))

    return max(combined_values)

### Pre-conditions: Takes a reference text, and integers k_mer and d...
### Post-conditions: Returns all most frequent k_mer patterns besides to their reverse complements with up to d mismatches in ref_text...
def Frequent_Words_With_Mismathces_And_Reverse_Complements(ref_text, k_mer, d):
    patterns = [] # will contain the most frequent k_mer patterns
    freq_map = {} # the frequency table of all k_mer patterns
    n = len(ref_text)

    # This for loop generates all k_mer patterns and their neighbors and populate the frequency table with them...
    for i in range(0, n-k_mer+1):
        pattern = ref_text[i:k_mer+i] 
        reverse_comlement_pattern = Reverse_Complement(pattern)
        neighborhood = Neighbors(pattern, d) # Generate all neighbors for that pattern as well it will contain the original pattern
        if n == k_mer:
            neighborhood = neighborhood.union(Neighbors(reverse_comlement_pattern, d))

        for neighbor in neighborhood:
            if freq_map.get(neighbor) == None :
                freq_map[neighbor] = 1
            else:
                freq_map[neighbor] += 1

    m = Max_Map_With_Reverse_Complements(freq_map) # Get the maximum frequency between all patterns
    # Find all patterns that have a maximum occurence in the text...
    for pattern in freq_map.keys():
        pattern_reverse_complement = Reverse_Complement(pattern)
        if freq_map[pattern] + freq_map.get(pattern_reverse_complement, 0) == m:
            patterns.append(pattern)
    
    return patterns
