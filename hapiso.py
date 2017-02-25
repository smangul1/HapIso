
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 04 18:45:18 2016

@author: Harry Yang
"""

import pysam
import sys
import numpy as np
import csv
from Bio import SeqIO
from scipy.cluster.vq import whiten


from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
import collections
from clustering import mean_shift_clustering
from collections import Counter
import scipy.stats as stat

debug_marker = False

if debug_marker is True:

    bam = '/Users/harryyang/Documents/Research/HapIso/set2primer0_rm_dup.bam'
    g1=35265595
    g2=35289548
    out='/Users/harryyang/Documents/Research/HapIso/set2primer2_summary.txt'
    genome = '/Users/harryyang/Documents/Research/HapIso/genome.fa'
    chr = "chr14"
    output = True
else:
    if len(sys.argv) < 6 :
        print "HapIso - Haplotype-specific Isoform reconstrunction"
        print "Written by Harry Yang and Serghei Mangul"
        print "[1] - The bam file to run HapIso"
        print "[2] - Left boundary of the Gene"
        print "[3] - Right boundary of the Gene"
        print "[4] - Output file"
        print "[5] - Chromosome"
        print "[6] - Reference Genome"
        print "If any question: email at harry2416@gmail.com"
        sys.exit(1)
    bam = sys.argv[1]
    g1 = int(sys.argv[2])
    g2 = int(sys.argv[3])
    out = sys.argv[4]
    chr = sys.argv[5]
    genome = sys.argv[6]
    output = True





"""
NOTE: The given gene coordinates do not matter: I supplied wrong coord but it finds it self so supplying the coordinate is unneccsary?
"""

"""
Result output init
"""
if output is True:
    result_name = bam.split('.bam')[0] + '_result.txt'
    result=file(result_name, 'w')
    orig_stdout = sys.__stdout__
    sys.stdout = result

SUMMARY_INFO = True

if SUMMARY_INFO == True:
    summary_info_file = out

def write_sumamry(string_to_write):
    with open(summary_info_file, 'ab') as f:
        f.write(string_to_write)

"""
Mapping
"""

samfile = pysam.Samfile(bam, "rb" )
print bam
print "Number of Mapped reads", samfile.mapped
print "Number of Unmapped reads", samfile.unmapped

# print g1, "----", g2

flag = 0

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=1)


readsA = []
for pileupcolumn in samfile.pileup("", g1, g2):
    if flag==0:
        g1=pileupcolumn.pos
        flag=1
    g2=pileupcolumn.pos
#print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
#sys.exit(1)


# print "New boundaries"
print g1,"----",g2

g1E=g1 # ENSEMBL
g2E=g2 # ENSEMBL


for pileupcolumn in samfile.pileup(chr, g1, g2):
    for pileupread in pileupcolumn.pileups:
        readsA.append(pileupread.alignment.qname)

readsA= np.unique(readsA)
# print readsA
# print "len(readsA)",len(readsA)


print "-------",chr


inFile = open(genome,'r')
for record in SeqIO.parse(inFile,'fasta'):
    # print record.id
    if record.id == chr:
        #print str(record.seq)[g1:g2]
        chr1Ref=str(record.seq)
        # print "l chr ",chr," ",len(chr1Ref)
        # print "------"
        break

n=g2-g1+1
m=len(array(readsA))





### BINARY MATRIX GENERATION

#data=np.array(range(n*m)).reshape((n, m))
#data[:]=-1
data=[]



read_list = []
for i in range(len(readsA)):
    read_list.append('N')

Matrix=[]



pos_list = []
for i in range(len(readsA)):
    pos_list.append(0)

base_call = []

# dataT=[]
# for i in range(len(readsA)):
#      dataT.append(0)
binary_matrix = []
marker = 0
for pileupcolumn in samfile.pileup(chr, g1, g2):
    dataT=[]
    for i in range(len(readsA)):
        dataT.append(0)
    if pileupcolumn.n > 1:
        ref_base = ""
        alt_base = ""
        ref_list = []
        alt_list = []
        for pileupread in pileupcolumn.pileups:

            itemindex = np.where(readsA==pileupread.alignment.qname)
            n=itemindex[0]
            if not pileupread.is_del and not pileupread.is_refskip: # to skip splicing and deletions ; tehre is still base but it comes from previous
                # print pileupread.query_position
                read_list[n]=pileupread.alignment.seq[pileupread.query_position]
                if pileupread.alignment.seq[pileupread.query_position].lower() == chr1Ref[pileupcolumn.pos].lower():
                    dataT[n]=1
                    ref_base = pileupread.alignment.seq[pileupread.query_position]
                    ref_list.append(ref_base)
                    #data[pileupcolumn.pos-g1][n]=0
                elif pileupread.alignment.seq[pileupread.query_position].lower() != chr1Ref[pileupcolumn.pos].lower():
                    dataT[n]=-1
                    alt_base = pileupread.alignment.seq[pileupread.query_position]
                    alt_list.append(alt_base)
                else:
                    print "NONE"
                    sys.exit(333)
                    #data[pileupcolumn.pos-g1][n]=1
            elif pileupread.query_position is None:
                dataT[n]=0
                # print "NONE "
                continue
            pos_list[n] = pileupread.query_position

        kN = 0
        k1 = 0
        k_n1 = 0
        #print dataT
        for i in range(len(readsA)):
            if dataT[i] == 0:
                kN = kN + 1
            elif dataT[i] == 1:
                k1 = k1 + 1
            elif dataT[i] == -1:
                k_n1 = k_n1 + 1
            else:
                exit(27)
        # print "N's=",kN
        # print "cov=",k1+k0
        # print "indel ",pileupcolumn.pos,"\t",pileupread.indel

        #print dataT
# """
# Cutoff parameters - currently 20% of the reads
# """

        if (k_n1 > len(readsA)*0.15 and k1 > 0.15*len(readsA)) or (k_n1 >= 0.90*(k1 + k_n1) and k_n1 > 0.30 * len(readsA)):
            data.append(dataT)
            alt_base_tuple = Counter(alt_list).most_common(1)
            ref_base_tuple = Counter(ref_list).most_common(1)
            # if reference is empty 
            if not ref_base_tuple:
                ref_base_tuple = [(chr1Ref[pileupcolumn.pos], 0)]           
            # if debug_marker:
            #     # print k_n1, k1
            #     # print ref_base_tuple
            #     print ref_list

            alt_base = alt_base_tuple[0]
            ref_base = ref_base_tuple[0]
            print ref_base, alt_base,
            base_pair = []
            base_pair.append(ref_base)
            base_pair.append(alt_base)
            base_call.append(base_pair)
            print pileupcolumn.pos
            #if myListPrev[n]!=pileupread.qpos & pileupcolumn.n>10:

            # dataT=[]
            # for i in range(len(readsA)):
            #     dataT.append(-1)
        binary_matrix.append(dataT)
        Matrix.append(read_list)
            #print myList
            #if myList.count('A')+myList.count('C')+myList.count('T')+myList.count('G')>10 :
                #print "@COV",pileupcolumn.pos+1,myList.count('A'),myList.count('C'),myList.count('T'),myList.count('G'),myList.count('N')
        read_list=[]
        for i in range(len(readsA)):
            read_list.append('N')

binary_matrix=array(binary_matrix).transpose()


"""
Variable Documentation
Base_call = n*2 matrix - n = # position that has alternative allele
row1 = ref_base
row2 = alt_base
"""



#samfile.close()

data1=array(data).transpose()





#print data[2035,:]
#print data[2036,:]
#print data[2037,:]

# print np.shape(data1)
# for i in range(len(data)):
#     print np.count_nonzero(data[i])=


"""
binary matrix checker
- finds the read with all non-zero entries
- If all the SNP positions are zero, it should be excluded from the clustering
"""


def binary_matrix_checker(matrix):
    """
    input - data - # col = # SNP positinos, # row = # reads
    """
    # print np.shape(matrix)
    # print matrix
    num_snp = len(matrix[0])
    num_read = len(matrix)
    new_matrix = []
    non_informative_matrix = []
    perfect_read_index = -1
    perfect_read_index_list = []
    # if debug_marker:
        # print "GFFDGDGQEWJFKLDLJFLKSDJFLKDSJLKDSJL"
        # print matrix
    for i in range(len(matrix)):
        # if debug_marker:
        #     if np.count_nonzero(matrix[i]) == num_snp:
        #         print "perfect_index",matrix[i], i 
        # for j in range(len(matrix[0])):
        if np.count_nonzero(matrix[i]) != 0:
            new_matrix.append(matrix[i])
            if np.count_nonzero(matrix[i]) == num_snp:
                if perfect_read_index == -1:
                    perfect_read_index = i
                perfect_read_index_list.append(matrix[i])
        else:
            non_informative_matrix.append(matrix[i])
    # most_freq = Counter(new_matrix.).most_common(1)
    # most_freq = stat.mode(new_matrix)
    # print most_freq
    
    # NEW
    new_perfect_read_list = []
    for item in perfect_read_index_list:
        # print item
        index_string = ",".join(str(x) for x in item)
        new_perfect_read_list.append(index_string)
    most_frequent = Counter(new_perfect_read_list).most_common(1)
    perfect_index = most_frequent[0][0]
    
    perfect_read_index = -1

    for i in range(len(matrix)):
        temp_index = ",".join(str(x) for x in matrix[i])
        # print temp_index
        if temp_index == perfect_index:
            perfect_read_index = i
            break

    # print "PERFECT INDEX ", perfect_read_index






    return new_matrix, non_informative_matrix, perfect_read_index



data1_pre_filter = data1
data1, data1_noinfo, read_for_haplo_index = binary_matrix_checker(data1)

# sys.exit(25)

# data1 = data1_pre_filter



# sys.exit(22)

"""
Clustering - Mean Shift
"""
# print "Clustering Started"
# if debug_marker:
#     print "DATA1" , data1
cluster_vector = mean_shift_clustering(data1)
# print "Clustering Ended", "\n Result : ", np.size(cluster_vector), np.count_nonzero(cluster_vector), cluster_vector
# print "Clustering Finished"
print "ASE: ", np.size(cluster_vector) - np.count_nonzero(cluster_vector), ":", np.count_nonzero(cluster_vector)
print "Base Calls", base_call

ASE = []
ASE.append(np.size(cluster_vector) - np.count_nonzero(cluster_vector))
ASE.append(np.count_nonzero(cluster_vector))
ase_name= bam.split('.bam')[0] + '_ASE.txt'
np.savetxt(ase_name, ASE)


"""
Haplotype call test # TEST
"""
haplo_one = []
haplo_two = []
# print len(data)
temp_data_hap_call=np.transpose(data1_pre_filter)
# print base_call
# print np.shape(temp_data_hap_call)
for i in range(len(temp_data_hap_call)):
    if temp_data_hap_call[i][read_for_haplo_index] == 1:
        haplo_one.append(base_call[i][0])
        haplo_two.append(base_call[i][1])
    elif base_call[i][0] == '':
        """
        if the SNP is shared by both haplotype
        """
        haplo_one.append(base_call[i][1])
        haplo_two.append(base_call[i][1])
    elif temp_data_hap_call[i][read_for_haplo_index] == -1:
        haplo_one.append(base_call[i][1])
        haplo_two.append(base_call[i][0])

    else:
        # print temp_data_hap_call
        # print data, i
        # print "EXIT 28"
        # sys.exit(28)
        continue

if debug_marker:
    print "Halpotype_marker:", data1_pre_filter[read_for_haplo_index]
    print "Haplotype_one is:", haplo_one, "Haplotype_two is", haplo_two
# print data1
# sys.exit(33)
# pca_matrix_output = open('./pca_matrix_output.txt','w')
binary_matrix_name = bam.split('.bam')[0] + '_binary_matrix.txt'
np.savetxt(binary_matrix_name, data1)

"""
Cluster Vector Recovery Step
"""
#
# for i in range()
# """
#  Visual Checking Step
# """
# for i in range(len(cluster_vector)):
#     print readsA[i].rsplit('/')[1:2], cluster_vector[i]
#






# sys.exit(23)





# print cluster_vector
# print "len(cluster_vector)", len(cluster_vector)

# to get correct ASE
if output:
    allele_one = 0
    allele_two = 0 
    for read in cluster_vector:
        if read == 0:
            allele_one += 1
        elif read == 1:
            allele_two += 1
        else:
            sys.exit(33333)
    print "ASE RATIO is : ", allele_one, ":", allele_two
    print "Haplotype_one is:", haplo_one, "Haplotype_two is", haplo_two
    sys.exit(33)


#list for allelese from cluster 0
a1=[]
a2=[]


"""
ERROR CORRECTION
"""
print pos_list
for pileupcolumn in samfile.pileup(chr, g1, g2):
    #print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
    #print chr1Ref[pileupcolumn.pos].lower()

    if pileupcolumn.n>=0:
        for pileupread in pileupcolumn.pileups:
            #if (pileupcolumn.pos>6694285):
            #print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
            #print '\tbase in read %s = %s' % (pileupread.alignment.qname, pileupread.alignment.seq[pileupread.qpos])

            #array.append(pileupread.alignment.seq[pileupread.qpos])

            #print pileupcolumn.pos,', '.join(array)
            #print "pileupread.alignment.qname ",pileupread.alignment.qname
            itemindex = np.where(readsA == pileupread.alignment.qname)
            #print "item=",itemindex[0]
            n = itemindex[0]
            #print n
            #print "PREVIOUS POSITION=", myListPrev[n]
            #print "current position=", pileupread.qpos
            # print n
            if pos_list[n] != pileupread.query_position:
                # to skip splicing and deletions ; there is still base but it comes from previous
                if n is None:
                    print "XXXXXX"
                if pileupread.query_position is None:
                    continue
                #print idx[n]
                if pileupread.alignment.seq[pileupread.query_position].lower() != 'n' and cluster_vector[n] == 0:
                    a1.append(pileupread.alignment.seq[pileupread.query_position].lower())
                elif pileupread.alignment.seq[pileupread.query_position].lower() != 'n' and cluster_vector[n] == 1:
                    a2.append(pileupread.alignment.seq[pileupread.query_position].lower())
                elif cluster_vector[n] == -10:
                    continue
                read_list = []
                for i in range(len(readsA)):
                    read_list.append('N')

                a1_a=a1.count('a')
                a1_c=a1.count('c')
                a1_t=a1.count('t')
                a1_g=a1.count('g')


                a2_a=a2.count('a')
                a2_c=a2.count('c')
                a2_t=a2.count('t')
                a2_g=a2.count('g')


                a1C='n'
                a1C_count=0
                a2C='n'
                a2C_count=0


def printme(a1_a, a1_c, a1_t, a1_g):
    "This prints a passed string into this function"
    print str
    return a1C

    if a1_a>a1_c and a1_a>a1_t and a1_a>a1_g:
        a1C='a'
        a1C_count=a1_a

    if a1_c>a1_a and a1_c>a1_t and a1_c>a1_g:
        a1C='c'
        a1C_count=a1_c

    if a1_t>a1_a and a1_t>a1_c and a1_t>a1_g:
        a1C='t'
        a1C_count=a1_t

    if a1_g>a1_a and a1_g>a1_c and a1_g>a1_t:
        a1C='g'
        a1C_count=a1_g


    #----------------
    if a2_a>a2_c and a2_a>a2_t and a2_a>a2_g:
        a2C='a'
        a2C_count=a2_a

    if a2_c>a2_a and a2_c>a2_t and a2_c>a2_g:
        a2C='c'
        a2C_count=a2_c

    if a2_t>a2_a and a2_t>a2_c and a2_t>a2_g:
        a2C='t'
        a2C_count=a2_t

    if a2_g>a2_a and a2_g>a2_c and a2_g>a2_t:
        a2C='g'
        a2C_count=a2_g


    #if a1C!=chr1Ref[pileupcolumn.pos].lower() or a1C!=a2C :
    if a1C!=a2C and a1C!='n' and a2C!='n':
        print pileupcolumn.pos,a1C,a1C_count,a2C,a2C_count,chr1Ref[pileupcolumn.pos].lower()
    a1=[]
    a2=[]







#data[idx==0][:,0]







MatrixC1=[]
MatrixC2=[]
ar=[]

print "len(Matrix[0])",len(Matrix[0])
print "len(Matrix)",len(Matrix)

for i in range(0,len(Matrix[0])):

    for j in range(0, len(Matrix)):
        ar.append(Matrix[j][i])

    if cluster_vector[i] == 1:
        MatrixC1.append(ar)
    else:
        MatrixC2.append(ar)
    ar=[]




print MatrixC1



#x = np.array([1,1,1,2,2,2,5,25,1,1])
#y = np.bincount(x)

print "MatrixC1 ", len(MatrixC1)
print "MatrixC2 ", len(MatrixC2)

print len(MatrixC1[0])

# sys.exit(1)


MatrixC1=array(MatrixC1)
MatrixC2=array(MatrixC2)


print len(MatrixC1[:,0])
print len(MatrixC1[0,:])


print len(MatrixC2[:,0])
print len(MatrixC2[0,:])

a=MatrixC1[:,0].tolist()
a[:] = (value for value in a if value != 'N')


import collections
counter=collections.Counter(a)
print(counter),"COUNTER"
# Counter({1: 4, 2: 4, 3: 2, 5: 2, 4: 1})
print(counter.values()),"COUNTER_VAL"
# [4, 4, 2, 1, 2]
print(counter.keys()),"COUNTER<KEYS"
# [1, 2, 3, 4, 5]
print(counter.most_common(1)),"COUNTER>MOSTCOMMON"
# [(1, 4), (2, 4), (3, 2)]


e='N'
for letter, count in counter.most_common(1):
    # print '%s: %7d' % (letter, count)
    e=letter

print "e=",e
if len(e) == 0:
    print "!!!!!"


haplo1=[]
for i in range(0,len(MatrixC1[0,:])):
    a=MatrixC1[:,i].tolist()
    a[:] = (value for value in a if value != 'N')
    counter=collections.Counter(a)
    counter.most_common(1)
    e='N'
    for letter, count in counter.most_common(1):
        e=letter
    # print i,letter
    haplo1.append(e)

print "haplo1"
print haplo1


a=MatrixC2[:,0].tolist()
a[:] = (value for value in a if value != 'N')


import collections
counter=collections.Counter(a)
print(counter)
# Counter({1: 4, 2: 4, 3: 2, 5: 2, 4: 1})
print(counter.values())
# [4, 4, 2, 1, 2]
print(counter.keys())
# [1, 2, 3, 4, 5]
print(counter.most_common(1))
# [(1, 4), (2, 4), (3, 2)]

e='N'
for letter, count in counter.most_common(1):
    # print '%s: %7d' % (letter, count)
    e=letter

print "e=",e
if len(e) == 0:
    print "!!!!!"



haplo2=[]
for i in range(0,len(MatrixC2[0,:])):
    a=MatrixC2[:,i].tolist()
    a[:] = (value for value in a if value != 'N')
    counter=collections.Counter(a)
    counter.most_common(1)
    e='N'
    for letter, count in counter.most_common(1):
        e=letter
    #print i,letter
    haplo2.append(e)


if len(haplo1) != len(haplo2):
    print "ERROR!"
    print "exit"
    # sys.exit(1)

print "len(haplo1)",len(haplo1)
print "len(haplo2)",len(haplo1)

if len(MatrixC1)>len(MatrixC2):
    ratio = float(len(MatrixC2))/len(MatrixC1)
else:
    ratio = len(MatrixC1)/float(len(MatrixC2))

print "ratio", ratio





for i in range(0, len(haplo1)):
    if i+g1+1>=g1E and i+g2+1<=g2E:
        if haplo1[i]!=haplo2[i] and haplo1[i]!='N' and haplo2[i]!='N':
            print "@SNP",i,i+g1+1,haplo1[i],haplo2[i], MatrixC1[:,i].tolist().count(haplo1[i]), MatrixC2[:,i].tolist().count(haplo2[i]),ratio
        elif haplo1[i]==haplo2[i] and haplo1[i]!='N' and haplo2[i]!='N':
            print "@HOMO",i,i+g1+1,haplo1[i],haplo2[i], MatrixC1[:,i].tolist().count(haplo1[i]), MatrixC2[:,i].tolist().count(haplo2[i]),ratio


print "MatrixC1 ", len(MatrixC1)
print "MatrixC2 ", len(MatrixC2)


print "ENSEMBL begin",g1E
print "ENSEMBL end ",g2E
print "ALLELE1", len(MatrixC1)
print "ALLELE2", len(MatrixC2)


haplo1E=[]
haplo2E=[]



#write haplotype staking into consideration coverage
for i in range(0, len(haplo1)):
    if i+g1+1>=g1E and i+g1+1<=g2E:
        # print i+g1+1,g1E,g2E
        # print MatrixC1[:,i].tolist().count(haplo1[i]), MatrixC2[:,i].tolist().count(haplo2[i])
        # print haplo1[i]
        # print haplo2[i]
        if  MatrixC1[:,i].tolist().count(haplo1[i])==1 and MatrixC2[:,i].tolist().count(haplo2[i])==1 :
            # print "both N"
            haplo1E.append('N')
            haplo2E.append('N')
        elif MatrixC1[:,i].tolist().count(haplo1[i])==1:
            # print "case1"
            haplo1E.append(haplo2[i])
            haplo2E.append(haplo2[i])
        elif MatrixC2[:,i].tolist().count(haplo2[i])==1:
            # print "case2"
            haplo2E.append(haplo1[i])
            haplo1E.append(haplo1[i])
        else:
            # print "Both good"
            haplo1E.append(haplo1[i])
            haplo2E.append(haplo2[i])
        if len(haplo1E)!=len(haplo2E):
            print haplo1E
            print haplo2E
            print "ERROR"
            # sys.exit(1)







print "@HAPLO1",''.join(haplo1E)
print "@HAPLO2",''.join(haplo2E)

print "len(H1)",len(haplo1)
print "len(H2)",len(haplo2)
print "len(H1)",len(haplo1E)
print "len(H2)",len(haplo2E)


print "done!"


def comparison(haplo_one, haplo_two):
    for i in range(len(haplo_one)):
        if haplo_one[i] != haplo_two[i] and haplo_one[i] != 'N' and haplo_two[i] != 'N':
            print "SNV"
            print "Haplo 1:", haplo_one[i]
            print "Haplo 2:", haplo_two[i]
            print i
        elif haplo_one[i] != haplo_two[i] and (haplo_one[i] == "N" or haplo_two[i] == 'N'):
            # print "INDEL"
            # print "Haplo 1:", haplo_one[i]
            # print "Haplo 2:", haplo_two[i]
            # print "POS:", i
            continue


comparison(''.join(haplo1E), ''.join(haplo2E))

"""
error correction
"""
def error_correction(haplo_one, haplo_two, cluster_vector, binary_matrix):
    if len(cluster_vector) != len(binary_matrix):
        sys.exit(34)
    corrected_reads = []
    for i in range(len(cluster_vector)):
        current_hap = []
        if cluster_vector[i] == 0:
            current_hap = haplo_one
        elif cluster_vector[i] == 1:
            current_hap = haplo_two
        corrected_read = ""

        for j in range(len(binary_matrix[i])):
            if binary_matrix[i][j] == 0:
                continue
            else:
                if current_hap[j] != 'N':
                    corrected_read = corrected_read + current_hap[j]
                else:
                    continue
        corrected_reads.append(corrected_read)
    return corrected_reads


corrected_reads = error_correction(''.join(haplo1E), ''.join(haplo2E), cluster_vector, binary_matrix)
corrected_list = []
for i in range(len(readsA)):
    dummy_line = readsA[i] + "  " + corrected_reads[i]
    corrected_list.append(dummy_line)

cor_name = bam.split('.bam')[0] + '_corrected_reads.txt'
np.savetxt(cor_name, corrected_list)









if output is True:
    sys.stdout = orig_stdout
    result.close()
