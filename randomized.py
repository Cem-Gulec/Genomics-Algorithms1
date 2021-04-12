import random
import timeit
from collections import Counter



# k: k-mer length
# generate sequence of length k
def generate_sequence(k, bases='ACGT'):
    return ''.join([random.choice(bases) for i in range(k)])

# n: how many to generate
# k: k-mer length
# this function mutates given consensus string
# then generates updated sequence lines
# where mutated string is placed
# example run: generate_mutated_sequence(10, 490)
def generate_mutated_sequence(n, k):   
    
    mutated_list = []
    
    # mutation locations for every 10 variation
    rand_mutate_locs = [random.sample(range(n), 4) for _ in range(n)]
    # random location to place mutated string in range [0,500)
    rand_place_locs  = [random.sample(range(k), 1) for _ in range(n)]

    with open('input.txt') as f: lines = f.readlines()

    # 10 different sequence line
    for i in range(n):
        consensus_list = list(consensus_str)
        line = list(lines[i])

        # applying mutation into consensus string
        for k in range(4):
            consensus_list[rand_mutate_locs[i][k]] = excluded_base_select(consensus_list[rand_mutate_locs[i][k]])
                
        # place mutated k-mer in a sequence line
        for j in range(10):
            line[int(rand_place_locs[i][0]) + j] = consensus_list[j]
        
        listToStr = ''.join([str(elem) for elem in line])
        mutated_list.append(listToStr)
    
    for x in mutated_list:
        print(x)

# among ACGT bases exclude selected base and pick random one
def excluded_base_select(exclude_base):
    rand_base = random.choice('ACGT')
    return excluded_base_select(exclude_base) if rand_base is exclude_base else rand_base 

def hamming_distance(seq1, seq2):
    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count = count + 1
    
    return count

def print_profile_list(profile_list):
    print("Profile matrix:")
    print("A:", end=" ")
    for i in range(len(profile_list)): print(profile_list[i][0], end=" | ")
    print("\nT:", end=" ")
    for i in range(len(profile_list)): print(profile_list[i][1], end=" | ")
    print("\nG:", end=" ")
    for i in range(len(profile_list)): print(profile_list[i][2], end=" | ")
    print("\nC:", end=" ")
    for i in range(len(profile_list)): print(profile_list[i][3], end=" | ")
    print("\n")

# k: length of k-mer
# for a sequence line return every variation possible
def search_all_variations(string, k):
    variations = []
    length = len(string)
    
    for i in range(length-k): variations.append(string[i:i+k])

    return variations

# returns 10 motifs with highest probability and their index
def select_new_motifs(prob_list):
    motif_list = []
    for line in prob_list:
        temp_list = []
        # max value in a line
        temp_list.append(max(line))
        # its index value
        possible_indexes = [i for i, e in enumerate(line) if e == max(line)]
        # append occurences of max value with different patterns
        for index in possible_indexes:
            temp_list.append(index)

        motif_list.append(temp_list)

    return motif_list        

# calculates probability values for every variation
def calculate_prob_variations(Dna, profile_list, k):    
    variation_list = []
    prob_list = []
    
    # obtain every possible variation from each sequence line
    # size: 10x491 (10 lines, 491 variation each)
    for line in Dna:
        variations_from_line = search_all_variations(line, k)
        variation_list.append(variations_from_line)
    
    # calculate probabilities according to profile matrix
    # for every variation
    for variation_line in variation_list:
        temp_list = []
        for variation in variation_line:
            prob_value = calculate_single_prob(profile_list, variation)
            temp_list.append(prob_value)
        prob_list.append(temp_list)
    
    return prob_list

# calculates given string's probability
# depending on the profile matrix values
def calculate_single_prob(profile_list, variation):
    # enumarating bases for later index usage
    base_dict = {"A":0, "T":1, "G":2, "C":3}
    
    # initially hold [0th col][base_column_index]
    prob_total = profile_list[0][base_dict[variation[0]]]
    # then continue on multiplicating
    for i in range(1, len(variation)):
        prob_total *= profile_list[i][base_dict[variation[i]]]
    
    return prob_total

def calculate_score(motif_list, consensus_str):    
    score = 0
    rows = len(motif_list)
    
    for i in range(rows):
        line = "".join(motif_list[i])
        score += hamming_distance(consensus_str, line)
    
    return score

def generate_consensus(motif_list):    
    k = len(motif_list[0])
    consensus_return = ""
    for i in range(k):
        column_list = []
        for j in range(10): 
            column_list.append(motif_list[j][i])
        base_counts = Counter(column_list)
        A_prob = base_counts['A']
        T_prob = base_counts['T']
        G_prob = base_counts['G']
        C_prob = base_counts['C']
        highest_val = max([A_prob, T_prob, G_prob, C_prob])
        index = [A_prob, T_prob, G_prob, C_prob].index(highest_val)
        
        if index   == 0: consensus_return += 'A'
        elif index == 1: consensus_return += 'T'
        elif index == 2: consensus_return += 'G'
        elif index == 3: consensus_return += 'C'
    
    return consensus_return

# k: motif length
def randomized_motif_search(k, Dna):
    
    start = timeit.default_timer()
    iteration_num = 1
    motif_list = []
    profile_list = []

    with open(Dna) as f: lines = f.readlines()

    # random locations to locate motifs
    rand_locs  = [random.sample(range(490), 1) for _ in range(10)]

    # select random motifs from each sequence line
    for i in range(10):
        line = list(lines[i])
        motif_list.append(line[rand_locs[i][0]:rand_locs[i][0]+k])

    # declaring consensus string from motif list
    consensus_str = generate_consensus(motif_list)
    print("consensus string: ", consensus_str)
    
    # calculating score for first matrix
    best_score = calculate_score(motif_list, consensus_str)  
    initial_score = best_score
    print("Iteration number:", iteration_num)
    print("consensus string: {}\n".format(consensus_str))
    print("Motif List:")
    print(*motif_list, sep="\n")
    print("\ninitial score: {}\n".format(best_score))

    while True:
        # generate profile matrix from selected motifs
        for i in range(k):
            column_list = []
            for j in range(10): column_list.append(motif_list[j][i])
            base_counts = Counter(column_list)
            A_prob = base_counts['A']/10
            T_prob = base_counts['T']/10
            G_prob = base_counts['G']/10
            C_prob = base_counts['C']/10
            profile_list.append([A_prob, T_prob, G_prob, C_prob])
    
        # calculating every possibilities probability for each sequence line
        prob_list = calculate_prob_variations(lines, profile_list, k)
    
        # obtaining new motifs for each sequence line
        # which are having highest probability values
        # P.S: if there are multiple variations with same highest probability value
        # then their combinations for motif matrix are tried
        # whichever grants the better score is selected among them
        highest_motif_values = select_new_motifs(prob_list)
        
        # calculated to see what is max number of occurence of same pattern
        max_length = 0
        for i in range(len(highest_motif_values)):
            if len(highest_motif_values[i]) > max_length:
                max_length = len(highest_motif_values[i]) - 1
        
        new_motif_list = []

        for i in range(10):
            index = highest_motif_values[i][1]
            new_motif_list.append(lines[i][index:index+k])

        # check if score reduced
        # if so replace previous motif matrix with these new ones
        # from these motifs calculate profil once again
        # until it does not reduce score anymore algorithm continues..
        if calculate_score(new_motif_list, consensus_str) < best_score:
            consensus_str = generate_consensus(new_motif_list)
            new_score = calculate_score(new_motif_list, consensus_str)
            print("****************************************************")
            print("Iteration number: ", iteration_num)
            print("\nNew motif list:", new_motif_list)
            print("\nnew score:\n", new_score)
            print_profile_list(profile_list)
            motif_list = new_motif_list[:]
            profile_list = []
            best_score = new_score
        else:
            print("initial score was:", initial_score)
            print("best score possible:", best_score)
            new_score = calculate_score(new_motif_list, consensus_str)
            print("stopped at score:", new_score)
            print("number of iterations passed:", iteration_num)
            break
        iteration_num += 1
    
    stop = timeit.default_timer()
    print("Runtime of Randomized Motif Search: {}".format(stop-start), "sec")
    return motif_list


if __name__ == "__main__":
    #  9-mer consensus string:   CGAGCATCC
    # 10-mer consensus string:  ACGAGCATCC
    # 11-mer consensus string: ACGAGCATCCT
    #consensus_str = "CGAGCATCC"
    
    # inputs: k-mer length, input file
    # output: best motif acquired according to algorithm
    randomized_motif_search(10, "input.txt")