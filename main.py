import random

def generate_string(N, alphabet='ACGT'):
    return ''.join([random.choice(alphabet) for i in range(N)])

if __name__ == "__main__":
    dna = generate_string(500)
    print(dna)