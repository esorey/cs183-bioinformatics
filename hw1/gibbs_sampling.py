import random
import numpy as np

def compute_pwm(kmers):
    '''Given a list of kmers, compute the pwm. Format is
    [{*char* : probability}]'''
    pwm = []
    k = len(kmers[0])
    n = len(kmers)
    for i in range(k):
        chars = [kmer[i] for kmer in kmers]
        prob_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for char in chars:
            prob_dict[char] += 1 / n
        pwm.append(prob_dict)
    return pwm

def pwm_score_on_kmer(pwm, kmer):
    '''Score a kmer against a pwm'''
    k = len(kmer)
    assert k == len(pwm)
    score = 1
    for j in range(k):
        char = kmer[j]
        score *= pwm[j][char]
    return score

def gibbs_sampling(seq_list, k, num_iters):
    '''Run the gibbs sampling algorithm for num_iters iterations on the given
    sequence list, using kmers. Return the pwm.'''
    # Select a random kmer from each of the sequences
    kmers = []
    for seq in seq_list:
        n_kmers = len(seq) - k + 1
        # Add a random kmer from this sequence to the list of kmers
        random_idx = random.choice(range(n_kmers))
        random_kmer = seq[random_idx:random_idx + k]
        kmers.append(random_kmer)


    for _ in range(num_iters):
        # Randomly leave out a kmer
        random_idx = random.choice(range(len(seq_list)))
        random_seq = seq_list[random_idx]
        leave_one_out_kmers = kmers[:random_idx] + kmers[random_idx+1:]

        # Compute the pwm from the remaining kmers
        pwm = compute_pwm(leave_one_out_kmers)

        # Compute the pwm score for each kmer in the left out sequence. Keep in
        # dictionary.
        kmer_scores = {}
        for i in range(len(random_seq) - k + 1):
            kmer = random_seq[i:i+k]
            if kmer in kmer_scores:
                kmer_scores[kmer] += pwm_score_on_kmer(pwm, kmer)
            else:
                kmer_scores[kmer] = pwm_score_on_kmer(pwm, kmer)

        # Normalize kmer scores to 1
        score_sum = sum(kmer_scores.values())
        for kmer in kmer_scores:
            kmer_scores[kmer] /= score_sum

        # Choose a kmer, weighted by their scores
        picked_kmer = np.random.choice(list(kmer_scores.keys()), 1, list(kmer_scores.values()))[0]

        # Put the selected kmer in the list, and recompute the pwm
        kmers = kmers[:random_idx] + [picked_kmer] + kmers[random_idx+1:]
        pwm = compute_pwm(kmers)

    return pwm

# Run the algorithm on the given sequences
with open('sequences.txt') as f:
    content = f.readlines()
content = [x.strip() for x in content]

# We'll use k=4
print(gibbs_sampling(content, 4, 10))
print(gibbs_sampling(content, 4, 100))
print(gibbs_sampling(content, 4, 1000))
