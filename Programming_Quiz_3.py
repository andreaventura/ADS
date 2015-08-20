
# coding: utf-8

import numpy
import time

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


# In[66]:

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


        


def readGenome(filename):
    genome=''
    with open(filename,'r') as fh:
        for line in fh:
            if not line[0]=='>':
                genome += line.rstrip()
    return genome
            


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



chr1=readGenome('chr1.GRCh38.excerpt.fasta')


import bisect

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

def test_extended_match(pmer, t,mm=2):
    mismatches=0
    for i in range(len(pmer)):
        if pmer[i]!= t[i]:
            mismatches+=1
            if mismatches > mm:
                return False, 4
    return True, mismatches


# In[72]:

def pigeonHoleEdit(p,genome,index, max_edit=4):
    k=index.k
    ival=index.ival
    span=index.span
    pmer=[]
    hits={}
    total_hits=0
    for i in range (4):
        pmer=p[i:i+span]
        index_hits=index.query(pmer)
        total_hits += len(index_hits)
        for hit in index_hits:
            start=hit-i #position of pmer hit
            if start >= 0 and start+len(p) <= len(genome): #again to make sure it is a legal start!
                edit_matrix=local_editDist(p, genome[start-max_edit: start+len(p)+max_edit])
                #get the edit distance and the position of the best match
                edits, idx = min((edits, idx) for (idx, edits) in enumerate(edit_matrix[-1]))
                if edits<=max_edit:
                    hits[start]=(idx,genome[start-max_edit:start+max_edit+len(p)], edits)
    return hits, total_hits


chr1=readGenome('chr1.GRCh38.excerpt.fasta')
index=SubseqIndex(chr1, 4, 2)


#question 1: What is the edit distance of the best match between  
# pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1?
p='GCTGATCGATCGTACG'
result=pigeonHoleEdit(p, chr1, index,3)[0]
editDistances=[]
for key in result:
    editDistances.append(result[key][2])
print 'the edit distance of the best match between pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1 is: {}'.format(min(editDistances))
    


# In[57]:

# Q2: What is the edit distance of the best match between 
# pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1?
p='GATTTACCAGATTGAG'
result=pigeonHoleEdit(p, chr1, index,3)[0]
editDistances=[]
for key in result:
    editDistances.append(result[key][2])
print 'the edit distance of the best match between pattern GATTTACCAGATTGAG and the excerpt of human chromosome 1 is: {}'.format(min(editDistances))

def official_overlap(a, b, min_length=3):
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

reads,qualities=readFastq("/Users/venturaa/Google Drive/Computation files/ADS_course/ERR266411_1.for_asm.fastq")


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


# In[ ]:

t0=time.time()
overlap=overlap_graph(reads,30)
t1=time.time()-t0
print overlap[0:3], t1

#Question 3: calculate number of reads overlaps (number of edges)
edges=[]
for node in overlap:
    edges.append(len(node))
print 'number of edges: {}'.format(sum(edges))


# Question 4: calculate the number of nodes with at least 1 edge
with_edges = 0
no_edges =0
for edge in edges:
    if edge > 0:
        with_edges += 1
    else: 
        no_edges += 1
        
print 'reads with edges:\t{} \nreads without edges:\t{}'.format(with_edges, no_edges)




