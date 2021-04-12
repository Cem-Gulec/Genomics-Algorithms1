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

# calculates probability values for every variation
def calculate_single_prob(profile_list, variation):
    # enumarating bases for later index usage
    base_dict = {"A": 0, "T": 1, "G": 2, "C": 3}

    # initially hold [0th col][base_column_index]
    prob_total = profile_list[0][base_dict[variation[0]]]
    # then continue on multiplicating
    for i in range(1, len(variation)):
        prob_total *= profile_list[i][base_dict[variation[i]]]

    return prob_total


def calculate_score(motif_list, consensus_str): #calculate total hamming distance between motiflist elements and consensus string
    score = 0
    rows = len(motif_list)

    for i in range(rows):
        score += hamming_distance(consensus_str, motif_list[i])

    return score

def roll(massDist): # roll a biased die for prob. list
    randRoll = random.uniform(0,sum(massDist))  #
    sumi = 0
    result = 1
    for mass in massDist:
        sumi += mass
        if randRoll < sumi:
            return result
        result += 1

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
def Gibbs_Sampler(k, input_file):
    start = timeit.default_timer()
    N=0 #total iteration number
    sequence_list = [] #list of dna sequence
    motif_list = [] # k-mer motif list
    profile_list = [] #profile table list
    motifs=[] # last updated motif list
    threshold=0 # value that break the iteration because of score of the algorithms no longer improve


    with open(input_file) as f: #read txt
        lines = f.readlines()

    for i in lines:
        i = i[:-1]
        sequence_list.append(list(i)) #txt to dna sequence list


    # random locations to locate motifs
    rand_locs = [random.sample(range(490), 1) for _ in range(10)]

    # select random motifs from each sequence line  #########STEP1#########
    for i in range(10):
        line = list(lines[i])
        motif_list.append(line[rand_locs[i][0]:rand_locs[i][0] + k])#########STEP1#########
    
    consensus_str = generate_consensus(motif_list)
    print("consensus string: ", consensus_str)
    best_score = calculate_score(motif_list, consensus_str)#best initial score before iterations
    print("best score initial:")
    print(best_score)


    while(N<1000): # 1000 iteration process on motif list

        #remove random motif from motif list
        random_item_from_list = random.choice(motif_list)
        removed_item_index = motif_list.index(random_item_from_list)
        motif_list.remove(random_item_from_list)

        # generate profile matrix from selected motifs

        for i in range(k):
            column_list = []
            for j in range(9): column_list.append(motif_list[j][i])
            base_counts = Counter(column_list)
            A_prob = base_counts['A']/ 9 # 9, because of remove one row from motif list
            T_prob = base_counts['T']/ 9
            G_prob = base_counts['G']/ 9
            C_prob = base_counts['C']/ 9
            profile_list.append([A_prob, T_prob, G_prob, C_prob])



        kmer_list=[]
        probablity_list=[] #probablities of motif list elements
        str1 = ""
        str1=str1.join(sequence_list[removed_item_index]) # full sequence line of removed motif from motif list

        while (i + k < len(str1) + 1): # take substring that have k length from sequence line one by one
            str2 = str1[i:i + k]
            i = i + 1
            kmer_list.append(str2)

        for i in kmer_list:

            probablity_list.append(calculate_single_prob(profile_list, i)) #calculate probablity for string



        indeks = roll(probablity_list)-1 # index come from die process

        for i in range(len(motif_list)):
            line= "".join(motif_list[i])
            motifs.append(line)

        motifs.insert(removed_item_index, kmer_list[indeks]) # add new motif to removed motif's position

        #print(calculate_score(motifs, consensus_str))


        motif_list = []
        for i in range(len(motifs)):
            line = list(motifs[i])
            motif_list.append(line)
        
        consensus_str = generate_consensus(motifs)
        print("consensus string: ", consensus_str)
        print("Final Motif List:",motifs)
        new_score = calculate_score(motifs, consensus_str)
        print("new score: ", new_score)
        if new_score < best_score:
            best_score = new_score
        else:
            threshold=threshold+1

        print("best score:",best_score)
        print(best_score)
        motifs = []
        N = N + 1

        if threshold==500: # stop when the best score repeats more then 500 times
            break
    stop = timeit.default_timer()
    print("Runtime of Randomized Motif Search: {}".format(stop-start), "sec")


if __name__ == "__main__":
    #  9-mer consensus string:   CGAGCATCC
    # 10-mer consensus string:  ACGAGCATCC
    # 11-mer consensus string: ACGAGCATCCT
    # consensus_str = "CGAGCATCC"  # 9luk 10luk ve 11lik için ayrı ayrı ayarlanıyor
    Gibbs_Sampler(9, "input.txt")