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

short_seqs = ['TTTATCCAATAATGGACACGTT', 'CATAAATTTCACAAAACATATG']

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

    # Use these exact matches as seeds to expand to near-matches on the entirety
    # of the sequence.
    match_regions = []
    for match_idx in seed_match_idxs:
        mismatch_score = 0.0
        pad_size = 1
        while mismatch_score < 0.5:
            candidate = long_seq[match_idx - pad_size:match_idx+4+pad_size]
            short_seq_substr = short_seq[pad - pad_size:pad+4+pad_size]
            if len(short_seq_substr) != len(candidate): # Our indices have wrapped
                break
            # Compute the new mismatch score
            mismatch_score = 0.0
            for i in range(len(candidate)):
                if candidate[i] != short_seq_substr[i]:
                    mismatch_score += 1 / len(candidate)
            pad_size += 1
        # Track the mismatch score, match candidate, and region
        match_regions.append((mismatch_score, candidate,
                            (match_idx - pad_size, match_idx + 4 + pad_size)))

    # Sort the possible matches by mismatch score
    sorted_regions = sorted(match_regions, key=lambda x: x[0])

    # Keep just the top 5
    best_regions = sorted_regions[:5]

    # Print the results
    print('Trying to match ' + short_seq + ' with the CASP1 sequence...')
    print('Results:')
    if match_regions == []:
        print('No matches!')
    else:
        for reg in best_regions:
            print("Match at region " + str(reg[2][0]) + ":" + str(reg[2][1]) +
                  " => " + long_seq[reg[2][0]:reg[2][1]])
    print('\n')

# Comments:
# Upon inspection, the best scoring matches (those printed first) appear to
# be the most likely matches.
