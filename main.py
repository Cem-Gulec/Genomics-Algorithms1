import random


# k: k-mer length
# generate sequence of length k
def generate_sequence(k, bases='ACGT'):
    return ''.join([random.choice(bases) for i in range(k)])

# n: how many to generate
# k: k-mer length
# this function mutates given consensus string
# then generates updated sequence lines
# where mutated string is placed
def generate_mutated_sequence(n, k):   
    
    consensus_str = "ACGAGCATCC"
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


if __name__ == "__main__":
    