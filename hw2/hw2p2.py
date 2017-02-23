def get_kmers(reads, k):
    '''Takes a list of reads and returns the list of kmers of those reads.'''
    # Break the seqs into kmers
    kmers = []
    for read in reads:
        n = len(read)
        for i in range(n - k + 1):
            kmers.append(read[i:i+k])
    return kmers

def assemble_kmers(kmers):
    '''Takes a list of kmers and assembles them using a de Bruijin graph and
    an eulerian walk.'''
    # The nodes of the de Bruijin graph are the distinct k-1mers, and the
    # edges go from left k-1mers to right k-1mers. We'll store the graph as
    # an adjacency list that also tracks the multiplicity of edges.
    nodes = []
    for kmer in kmers:
        left = kmer[:-1]
        right = kmer[1:]
        nodes.append(left)
        nodes.append(right)
    nodes = list(set(nodes)) # hack to remove duplicates

    # Build the adjacency list.
    adj_list = {n : [] for n in nodes}
    for kmer in kmers:
        left = kmer[:-1]
        right = kmer[1:]
        adj_list[left].append(right)

    # Get the eulerian path
    assembled = assemble_de_bruijin(adj_list)
    return assembled

def assemble_de_bruijin(adj_list):
    '''Given an adjacency list of a de bruijin graph, find an eulerian path
    through the graph. Return the assembled genome.'''
    # First find the starting point
    seen = []
    for key in adj_list:
        seen.extend(adj_list[key])
    start = list((set(adj_list.keys()) - set(seen)))[0]
    res = start
    curr = start
    num_edges = sum(map(len, adj_list.values()))
    for i in range(num_edges):
        if adj_list[curr] == []:
            break
        res += adj_list[curr][0][-1]
        curr = adj_list[curr][0]
    return res

# Use the functions to do the assembly
with open('PS2_Q2_Sequence') as inp:
    reads = inp.readlines()
reads = [x.strip() for x in reads] # Remove newlines

thirty_mers = get_kmers(reads, 30)
assembled_string = assemble_kmers(thirty_mers)
print("The assembled string is: ")
print(assembled_string)
print("Its length is: " + str(len(assembled_string)))
