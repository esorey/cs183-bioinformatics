class SuffixArray:
    '''A simple implementation of suffix arrays.'''

    def __init__(self, string):
        '''Take a string and generate the array.'''
        self.string = string + '$' # Add a sentinel character to the end of the string.
        self.init_array()

    def init_array(self):
        '''Build the suffix array.'''
        # Build a dictionary of the suffixes and their offsets.
        suffix_dict = {}
        for i in range(len(self.string)):
            suffix_dict[self.string[i:]] = i

        # Sort the suffixes
        sorted_suffixes = sorted(suffix_dict.keys())

        # Put the offsets in the suffix array
        suffix_arr = []
        for i in range(len(sorted_suffixes)):
            suffix_arr.append((sorted_suffixes[i], suffix_dict[sorted_suffixes[i]]))
        self.arr = suffix_arr

    # def query(self, query):
    #     '''Query the suffix array for a given string. Return a list of indices corresponding
    #        to all locations of matches with the query.'''
    #     n = len(self.arr)
    #     mid_idx = 
    #
    #
    # def bin_search_leftmost(self, query):
    #     '''Given a suffix array of the form [('bar', 3), ('foo', 2), ...], this returns the
    #     index of the leftmost entry whose suffix contains the given query.'''
