def get_kmer_freqs_and_nframes(seq_list, k):
    kmers_dict = {}
    total_kmer_frames = 0
    for seq in seq_list:
        for i in range(len(seq) - k + 1):
            total_kmer_frames += 1
            kmer = seq[i:i+k]
            if kmer in kmers_dict:
                kmers_dict[kmer] += 1
            else:
                kmers_dict[kmer] = 1
    kmers_dict = {key: kmers_dict[key] for key in kmers_dict}
    return kmers_dict, total_kmer_frames

def get_common_kmers(seq_list, k):
    kmers_dict, nframes = get_kmer_freqs_and_nframes(seq_list, k)
    res = []
    expected_freq = nframes * (1 / 4) ** k
    print(expected_freq)
    for kmer in kmers_dict:
        if kmers_dict[kmer] >= expected_freq:
            res.append(kmer)
    return res
