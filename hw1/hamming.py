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

def hamming_distance_same_size(seq1, seq2):
    '''Compute the hamming distance between the two given sequences of the same
    length.'''
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist

def min_hamming_distance(seq, pattern):
    '''Compute the minimum hamming distance between pattern and any subsequence
    of seq.'''
    assert len(seq) >= len(pattern)
    k = len(pattern)
    min_dist = len(pattern) + 1 # initialize
    min_dist_subseq = ''
    for i in range(len(seq) - k + 1):
        subseq = seq[i:i+k]
        dist = hamming_distance_same_size(pattern, subseq)
        if dist < min_dist:
            min_dist = dist
            min_dist_subseq = subseq
    return (min_dist_subseq, min_dist)

def min_hamming_distances_seqs(seq_list, pattern):
    ''' Apply min_hamming_distance to each sequence in seq_list; return the
    results as a list.'''
    res = []
    for seq in seq_list:
        res.append(min_hamming_distance(seq, pattern))
    return res


with open('sequences.txt') as f:
    content = f.readlines()
content = [x.strip() for x in content]

patterns = ['TTGTAGG', 'GAGGACC', 'TATACGG', 'CCGCAGG', 'CAGCAGG']

score_dict = {p: None for p in patterns}
sub_seq_dict = {p: None for p in patterns}
for p in patterns:
    score_dict[p] = sum([x[1] for x in min_hamming_distances_seqs(content, p)])
    sub_seq_dict[p] = [x[0] for x in min_hamming_distances_seqs(content, p)]

# The most likely motif is 'CAGCAGG'. Let's find its entropy
motif = 'CAGCAGG'
subseqs = sub_seq_dict[motif]
pwm = compute_pwm(subseqs)
entropy = 0
for i in range(len(motif)):
    entropy -= pwm[i][motif[i]] * np.log(pwm[i][motif[i]]) # use natural log
print(entropy)
