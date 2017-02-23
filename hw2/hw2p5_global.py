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


# Initialize scores for global alignment
backtrack_table[0][0] = 'start'
for i in range(1, len(long_seq)+1):
    dp_table[0][i] = -1 * i
    backtrack_table[0][i] = 'l'
for i in range(1, len(short_seq)+1):
    dp_table[i][0] = -1 * i
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
        max_score = max([skip_long_score, skip_short_score, align_score])

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
        backtrack_table[i][j] = backtrack


# Reconstruct the optimal alignment
optimal_alignment_score = dp_table[-1][-1] # Very last element
alignment_steps = ''
i = len(short_seq)
j = len(long_seq)
while backtrack_table[i][j] != 'start':
    curr = backtrack_table[i][j]
    alignment_steps += curr[0]
    if curr[0] == 'l':
        j -= 1
    elif curr[0] == 'u':
        i -= 1
    elif curr[0] == 'd':
        j -= 1
        i -= 1

# Reverse the alignment steps
alignment_steps = alignment_steps[::-1]

# Step through the alignment steps and build up the alignment. Since we're
# going from top to bottom, all of the backtrack entries must be reversed (ie
# left becomes right, etc)
long_alignment = ''
short_alignment = ''

i_short = 0
i_long = 0
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

# Print the alignment to a file
with open('p5_global_output.txt', 'w') as outfile:
    outfile.write(long_alignment)
    outfile.write('\n')
    outfile.write(short_alignment)

# Comments:
# Looking at the file, we see that the global alignment works by inserting
# skips into the longer sequence so that it exactly matches the shorter
# sequence.
