''' A series of functions I created while studying the Advances DNA Sequencing
course on Coursera
'''
import bisect
import numpy
import time
import bm_preproc

def reverseComplement(t):
    translation={'A':'T', 'T':'A', 'G': 'C', 'C':'G'}
    revcom=''
    for i in range(len(t)):
        revcom=translation[t[i]]+revcom
    return revcom


def readFastq(filename):
    sequences=[]
    qualities=[]
    with open(filename,'r') as fh:
        while True:
            fh.readline()
            seq=fh.readline().rstrip()
            fh.readline()
            qual=fh.readline().rstrip()
            if len(seq)==0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def readGenome(filename):
    genome=''
    with open(filename,'r') as fh:
        for line in fh:
            if not line[0]=='>':
                genome += line.rstrip()
    return genome


def naive(p,t):
    hits=[]
    rev_p=reverseComplement(p)
    for i in range(len(t)-len(p)+1):
        match_F=True
        match_R=True
        #check forward
        for j in range (len(p)):
            if not t[i+j]==p[j]:
                match_F=False
                break
        #check reverse
        for j in range (len(p)):
            if not t[i+j]==rev_p[j]:
                match_R=False
                break
        if match_F or match_R:
            hits.append(i)
    return hits

def boyer_moore(p,p_bm,t):
    """Do Boyer-Moore matching. p=pattern, t=text, p_bm=BoyerMoore object for p"""
    i=0
    occurrences=[]
    while i < len (t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range (len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(skip_bc, skip_gs, shift)
                mismatched=True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs=p_bm.match_skip()
            shift=max(shift,skip_gs)
        i += shift
    return occurrences

class Index(object):
    def __init__(self, t, k):
        self.k=k
        self.index=[]
        for i in range(len(t)-k+1):
            self.index.append((t[i:i+k], i))
        self.index.sort()

    def query(self, p):
        kmer=p[:self.k]
        #print kmer
        i=bisect.bisect_left(self.index, (kmer, -1))
        hits=[]
        #print i
        while i < len(self.index):
            if self.index[i][0] != kmer:

                break
            hits.append(self.index[i][1]) #appends the position of the kmer hit in the index
            i += 1
        return hits

class SubseqIndex(object):
    """Holds a subsequence index for a text T"""

    def __init__(self, t, k, ival):
        self.k=k
        self.ival=ival
        self.index=[]
        self.span=(k-1)*ival+1 #number of characters spanned by a kmer given ival
        for i in range (len(t)-self.span+1):
            self.index.append((t[i:i+self.span:ival], i))
        self.index.sort()

    def query(self, p):
        subseq=p[:self.span:self.ival]
        #print subseq
        i=bisect.bisect_left(self.index, (subseq, -1))
        hits=[]
        while i < len(self.index):
            if self.index[i][0] != subseq:

                break
            hits.append(self.index[i][1]) #appends the position of the kmer hit in the index
            i += 1
        return hits

def pigeon_hole(p, index, genome, mm=2):
    pmer=[] #creates list n of kmers from the pattern
    n=3
    k=8
    results={}
    total_hits=0
    for i in range (n):
        pmer=p[i*k : i*k+k]
        #print pmer
    #now check if there are matches in the indext
        hits=index.query(pmer)
        total_hits+=len(hits)
        #extend if there are hits
        if len(hits)>0:
            for hit in hits:
                start=hit-(i*k) #calculates the predicted start of the match
                match=True
                if start >= 0 and start not in results.keys(): #makes sure the start is legal and a match at that start has not been found yet
                    mismatches=0
                    for j in range(len(p)):
                        if p[j] != genome[start+j]:
                            mismatches+=1
                        if mismatches > mm:
                            match=False
                            break
                    if match:
                        results[start]=(mismatches, genome[start: start+len(p)])
    return results,total_hits

def test_extended_match(pmer, t,mm=2):
    mismatches=0
    for i in range(len(pmer)):
        if pmer[i]!= t[i]:
            mismatches+=1
            if mismatches > mm:
                return False, 4
    return True, mismatches

def pigeon_hole_sub(p,genome,index, mm=2):
    k=index.k
    ival=index.ival
    span=index.span
    print 'span', span, k, ival
    pmer=[]
    hits={}
    total_hits=0
    for i in range (3):
        pmer=p[i:i+span]
        #print pmer
        index_hits=index.query(pmer)
        total_hits += len(index_hits)
        for hit in index_hits:
            #print hit
            start=hit-i
            if start >= 0 and start+span <= len(genome): #again to make sure it is a legal start!
                match, mismatches=test_extended_match(pmer, genome[start: start+len(pmer)], mm)
                if match:
                    hits[start]=(genome[start:start+len(pmer)], mismatches)
    return hits, total_hits

def getDelta(x,y):
    return x!=y
    #if x==y:
    #    return 0
    #return 1

def globalEditDist(a,b):
    start_time=time.time()
    mat=numpy.zeros((len(a)+1, len(b)+1), int)
    #initialize matrix
    for i in range (len(a)+1):
        mat[i,0]=i
    for i in range (len(b)+1):
        mat[0,i]=i
    for row in range(1,len(a)+1):
        for col in range(1,len(b)+1):
            delt=getDelta(a[row-1], b[col-1])+mat[row-1,col-1]
            ins_a=mat[row,col-1]+1
            ins_b=mat[row-1,col]+1
            #print 'delta=%d insert_row=%d, insert_col=%d' %(delt, ins_a, ins_b)
            #now insert minimum valuea
            mat[row,col]=min(delt,ins_a,ins_b)
    #format the array so that it also has the

    return mat


def local_editDist(a,b):
    start_time=time.time()
    mat=numpy.zeros((len(a)+1, len(b)+1), int)
    #initialize matrix
    for i in range (len(a)+1):
        mat[i,0]=i
    for i in range (len(b)+1):
        mat[0,i]=0


    for row in range(1,len(a)+1):
        for col in range(1,len(b)+1):
            delt=getDelta(a[row-1], b[col-1])+mat[row-1,col-1]
            ins_a=mat[row,col-1]+1
            ins_b=mat[row-1,col]+1
            #print 'delta=%d insert_row=%d, insert_col=%d' %(delt, ins_a, ins_b)
            #now insert minimum valuea
            mat[row,col]=min(delt,ins_a,ins_b)
    #format the array so that it also has the

    return mat

def pigeonHoleEdit(p,genome,index, max_edit=4):
    k=index.k
    ival=index.ival
    span=index.span
    #print 'span', span, k, ival
    pmer=[]
    hits={}
    total_hits=0
    for i in range (4):
        pmer=p[i:i+span]
        #print pmer
        index_hits=index.query(pmer)
        total_hits += len(index_hits)
        #print total_hits
        for hit in index_hits:
            #print hit
            start=hit-i #position of pmer hit
            if start >= 0 and start+len(p) <= len(genome): #again to make sure it is a legal start!
                edit_matrix=local_editDist(p, genome[start-max_edit: start+len(p)+max_edit])
                #get the edit distance and the position of the best match
                edits, idx = min((edits, idx) for (idx, edits) in enumerate(edit_matrix[-1]))
                if edits<=max_edit:
                    hits[start]=(idx,genome[start-max_edit:start+max_edit+len(p)], edits)
    return hits, total_hits

def official_overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    # if len(b)< min_length:
    #    return 0
    #if len(a) < min_length:
    #    return 0
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def indexSuffix(reads, k=30):
    suffx_dic={}
    for read in reads:
        suffix=read[-k:] #get the suffix of length k
        #add the suffix as a key to the dictionary
        suffx_dic[suffix]=set()
        #populate with the corresponding kmers
    for i in range(len(reads)):
        read = reads[i]
        for j in range (len(read)-k+1): #remove '+1' if you want to avoid to add the suffix itself
            kmer=read[j:j+k]
            if kmer in suffx_dic:
                suffx_dic[kmer].add(i)
    return suffx_dic

def overlap_graph(reads, k=30):
    dic = indexSuffix(reads, k)  # first create the suffix dictionary and populate.
    ''' each key in dic is the suffix of at least one read. to each key is associated a set of indexes
    pointing to the reads in which a kmer is identical to the key. this way, given a suffix it is easy to know in which
    reads to look for an overlap
    and which reads can be ignored.'''

    #create a list in which each index corresponds to a read, and the associated set corresponds to the overlapping reads
    graph = []
    read_index=0
    for read in reads:
        overlaps=[]
        suffix=read[-k:]
        candidate_reads=dic[suffix] #indexes of the reads that contain the sequence corresponding to the suffix
        for candidate in candidate_reads:
            if read_index != candidate:   #to avoid finding overlaps to itself!
                hit=official_overlap(read, reads[candidate], k)
                if hit > 0:
                    overlaps.append((candidate, hit)) #tuple containing the index of the read and the length of the match
                    #print 'overlap found between: \nread %d %s\nread %d %s' %(read_index,reads[read_index],candidate,reads[candidate])
        graph.append(overlaps) #add to the graph the overlaps found (if any) for the read examined
        read_index += 1
    return graph
