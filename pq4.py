from adsFunctions import *

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=20)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest


def scs_all(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        #print sup
        if len(shortest_sup) > 0: #prevents this test from running when there is no shortest sup yet
            if len(sup) == len (shortest_sup[0]): #if the superstring analyzed is exactly as long as the longest so far add it to result list and
                shortest_sup.append(sup)
        if len(shortest_sup)==0 or len(sup) < len(shortest_sup[0]): #found new shortest superstring
            del shortest_sup[:]
            shortest_sup=[sup]  #empty the list and assign sup as the new first element

    return shortest_sup  # return shortest

'''Question 4: How many As are there in the full, assembled genome?'''
# let's try with the brute force approach first

def greedy_scs(ss, min_len=40):
    ''' this function takes an array of strings and tries to find the shortest common superstring using the greedy approach '''
    graph=overlap_graph(ss, min_len)
    # on approach is to progressively collapse the reads with the longest overlaps and then repeat
    while True:
        overlaps=[]
        for node in graph:
            if len(node)>0:
                max_overlap=numpy.argmax([olap[1] for olap in node])#the sequence with greatest overlap with each node
                overlaps.append(max_overlap)
            else: overlaps.append([])

        go_index=numpy.argmax(overlaps) #get the index of the node with the greatest overlap
        node=graph[go_index] #node having the greatest overlap
        go_seq_a=ss[go_index] #5'seq
        go_seq_b_index=numpy.argmax([for x[1] in x in node])
        merged_seq=ss[go_index]




    return ograph

    #test
