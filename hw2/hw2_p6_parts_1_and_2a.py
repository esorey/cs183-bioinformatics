# Part 1
# The alignment score is 12.

# Part 2a
with open('PS2_Q5_Sequence.txt') as infile:
    long_seq = infile.readlines()
long_seq = [x.strip() for x in long_seq[1:]]
temp = ''
for s in long_seq:
    temp += s
long_seq = temp

short_seqs = ['TCAGGTCACTCCATGCACAT', 'CAGTTCTGATTCTTTAATGG', 'AACTCAAG',
              'CATTAATT']

for short_seq in short_seqs:
    n = len(short_seq)
    pad = int((n - 4) / 2)
    # Get the middle 4mer
    middle4 = short_seq[pad:pad+4]

    # Find exact match indices in the long sequence.
    seed_match_idxs = []
    # Iterate over all 4mers
    for i in range(len(long_seq) - 3):
        substring = long_seq[i: i + 4]
        if substring == middle4:
            seed_match_idxs.append(i)

    # Use these exact matches as seeds to identify exact matches on the entirety
    # of the short sequence
    match_regions = []
    for match_idx in seed_match_idxs:
        candidate = long_seq[match_idx-pad:match_idx+pad+4]
        if candidate == short_seq:
            match_regions.append((match_idx - pad, match_idx + pad + 4))

    # Print the results
    print('Trying to match ' + short_seq + ' with the CASP1 sequence...')
    print('Results:')
    if match_regions == []:
        print('No matches!')
    else:
        for reg in match_regions:
            print("Match at region " + str(reg[0]) + ":" + str(reg[1]) + "!")
    print('\n')
