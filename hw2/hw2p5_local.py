with open('PS2_Q5_Sequence.txt') as infile:
    long_seq = infile.readlines()
long_seq = [x.strip() for x in long_seq[1:]]
temp = ''
for s in long_seq:
    temp += s
long_seq = temp

short_seq = 'ATCTCAAACACATGCGGGACCCCAGATA'

dp_table = []
for i in range(len(short_seq) + 1):
    dp_table.append([0] * (len(long_seq) + 1))
# backtrack encoding: 'l' -> came from left, 'u' -> came from up, 'd' ->
# came from diagonal, 'r' -> reset (local alignment only)
backtrack_table = []
for i in range(len(short_seq) + 1):
    backtrack_table.append([0] * (len(long_seq) + 1))


# Initialize scores for local alignment
backtrack_table[0][0] = 'start'
for i in range(1, len(long_seq)+1):
    dp_table[0][i] = 0
    backtrack_table[0][i] = 'l'
for i in range(1, len(short_seq)+1):
    dp_table[i][0] = 0
    backtrack_table[i][0] = 'u'

# Fill out the dynamic programming table, going across rows and then down
# columns
MISMATCH = -1
MATCH = 1
SKIP = -1

for i in range(1, len(short_seq)+1):
    for j in range(1, len(long_seq)+1):
        skip_long_score = SKIP + dp_table[i][j-1]
        skip_short_score = SKIP + dp_table[i-1][j]
        if long_seq[j-1] == short_seq[i-1]:
            align_score = MATCH + dp_table[i-1][j-1]
        else:
            align_score = MISMATCH + dp_table[i-1][j-1]
        max_score = max([0, skip_long_score, skip_short_score, align_score])

        # Update the score
        dp_table[i][j] = max_score

        # Update the backtracking
        backtrack = ''
        if max_score == skip_long_score:
            backtrack += 'l'
        if max_score == skip_short_score:
            backtrack += 'u'
        if max_score == align_score:
            backtrack += 'd'
        if max_score == 0:
            backtrack += 'r'
        backtrack_table[i][j] = backtrack


# Reconstruct the optimal alignment
# First, find the highest score out of the dp_table, along with the indices
# at which it occurs.
max_score = max([max(row) for row in dp_table])
max_idxs = []
for i in range(len(short_seq)+1):
    for j in range(len(long_seq)+1):
        if dp_table[i][j] == max_score:
            max_idxs.append((i,j))

# just use the first one; for these seqs, the max only occurs once anyway
max_idx = max_idxs[0]

# Backtrack from this max score until we find a score of 0 (ie, the start
# point for the local alignment)
i = max_idx[0]
j = max_idx[1]
alignment_steps = ''
while dp_table[i][j] != 0:
    step = backtrack_table[i][j][0]
    alignment_steps += step
    if step == 'l':
        j -= 1
    elif step == 'u':
        i -= 1
    elif step == 'd':
        i -= 1
        j -= 1


# Reverse the alignment steps
alignment_steps = alignment_steps[::-1]

# Step through the alignment steps and build up the alignment. Since we're
# going from top to bottom, all of the backtrack entries must be reversed (ie
# left becomes right, etc)
long_alignment = '-' * j # skip to where the alignment begins
short_alignment = '-' * i # skip to where the alignment begins

i_short = i
i_long = j
for step in alignment_steps:
    if step == 'l':
        short_alignment += '-'
        i_long += 1
    elif step == 'u':
        long_alignment += '-'
        i_short += 1
    elif step == 'd':
        long_alignment += long_seq[i_long]
        short_alignment += short_seq[i_short]
        i_long += 1
        i_short += 1
# Pad out with skips ('-')
while i_short < len(short_seq):
    short_alignment += '-'
    i_short += 1
while i_long < len(long_seq):
    long_alignment += '-'
    i_long += 1

# Print the alignment to a file
with open('p5_local_output.txt', 'w') as outfile:
    outfile.write(long_alignment)
    outfile.write('\n')
    outfile.write(short_alignment)

# Comments:
# This local alignment results in finding a substring in the long sequence
# ('ATCTCAAACACTTGCGTGACCTCAGATA') that is very close to exactly matching
# the entirety of the short sequence ('ATCTCAAACACATGCGGGACCCCAGATA'). This
# differs from the global alignment in that it does not use skips; rather it
# finds (in this case) the best matching substring within the larger sequence.
