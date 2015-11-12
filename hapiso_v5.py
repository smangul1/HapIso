## This script contains some example code 
## illustrating ways to to use the pysam 
## interface to samtools.
##
## The unit tests in the script pysam_test.py
## contain more examples.
##

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


if len(sys.argv)<4 :
    print "[1]-bam"
    print "[2]-left boundary"
    print "[3]-right boundary"
    print "[4] -out file"
    print "[5] - chr"
    sys.exit(1)





# data generation
#data = vstack((rand(150,2) + array([.5,.5]),rand(150,2)))

# computing K-Means with K = 2 (2 clusters)
#centroids,_ = kmeans(data,2)
# assign each sample to a cluster
#idx,_ = vq(data,centroids)





data3=array([[ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.]])

data2=[]

data2.append([ 0.,  0.,  0.,  0.,  0.])
data2.append([ 0.,  2.,  1.,  0.,  1.])
data2.append([ 0.,  2.,  1.,  0.,  1.])



data4=array(data2)

data7=np.transpose(data4)



print data7
print data7[1,:]
print max(set(data7[1,:].tolist()), key=data7[1,:].tolist().count)



# data generation
#data = vstack((rand(5,2) + array([.5,.5]),rand(5,2)))

print data4
print "----"
# computing K-Means with K = 2 (2 clusters)
centroids,_ = kmeans(data4,2)
# assign each sample to a cluster
idx,_ = vq(data4,centroids)

print "----"
print idx
print "----"



bam=sys.argv[1]
g1=int(sys.argv[2])
g2=int(sys.argv[3])
out=sys.argv[4]
chr=sys.argv[5]


samfile = pysam.Samfile(bam, "rb" )
print "mapped=",samfile.mapped
print "UNmapped=",samfile.unmapped

print g1,"----",g2

flag=0;

readsA = []
for pileupcolumn in samfile.pileup("", g1, g2):
    if flag==0:
        g1=pileupcolumn.pos
        flag=1
    g2=pileupcolumn.pos
#print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
#sys.exit(1)


print "New boundaries"
print g1,"----",g2

g1E=g1 # ENSEMBL
g2E=g2 # ENSEMBL


for pileupcolumn in samfile.pileup(chr, g1, g2):
    for pileupread in pileupcolumn.pileups:
        readsA.append(pileupread.alignment.qname)

readsA= np.unique(readsA)
print readsA
print "len(readsA)",len(readsA)


print "-------",chr


inFile = open('/u/home/s/serghei/project/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa','r')
for record in SeqIO.parse(inFile,'fasta'):
    print record.id
    if record.id == chr:
        #print str(record.seq)[g1:g2]
        chr1Ref=str(record.seq)
        print "l chr ",chr," ",len(chr1Ref)
        print "------"
        break

n=g2-g1+1
m=len(array(readsA))


#data=np.array(range(n*m)).reshape((n, m))
#data[:]=-1
data=[]



myList=[]
for i in range(len(readsA)):
    myList.append('N')

Matrix=[]


myListPrev=[]
for i in range(len(readsA)):
    myListPrev.append(0)



dataT=[]
for i in range(len(readsA)):
     dataT.append(-1)


for pileupcolumn in samfile.pileup(chr, g1, g2):
    print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
    
    print chr1Ref[pileupcolumn.pos].lower()
    
    if pileupcolumn.n>=0:
        
        for pileupread in pileupcolumn.pileups:
            #print pileupread
            #if (pileupcolumn.pos>6694285):
            #print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
                
    
            #print pileupcolumn.pos,', '.join(array)
            #print "pileupread.alignment.qname ",pileupread.alignment.qname
            itemindex = np.where(readsA==pileupread.alignment.qname)
            #print "item=",itemindex[0]
            n=itemindex[0]
            #print n
            #print "PREVIOUS POSITION=", myListPrev[n]
            print "current position=", pileupread.query_position
            
            if not pileupread.is_del and not pileupread.is_refskip and pileupread.query_position!='None': # to skip splicing and deletions ; tehre is still base but it comes from previous
                print pileupread.query_position
                myList[n]=pileupread.alignment.seq[pileupread.query_position]
                if pileupread.alignment.seq[pileupread.query_position].lower() == chr1Ref[pileupcolumn.pos].lower() :
                    dataT[n]=0
                    #data[pileupcolumn.pos-g1][n]=0
                else:
                    dataT[n]=1
                    #data[pileupcolumn.pos-g1][n]=1
            myListPrev[n]=pileupread.query_position


        kN=0
        k1=0
        k0=0
        #print dataT
        for i in range(len(readsA)):
            if dataT[i]==-1:
                kN=kN+1
            elif dataT[i]==1:
                k1=k1+1
            else:
                k0=k0+1
    #print "N's=",kN
    #print "cov=",k1+k0
    #print "indel ",pileupcolumn.pos,"\t",pileupread.indel
        
        #print dataT
    
        if k0>5 and k1>1:
            for i in range(len(readsA)):
                if dataT[i]==-1:
                    dataT[i]=0
            print pileupcolumn.pos,"    ",kN," ",k0,"  ",k1
            #print dataT
            data.append(dataT)
                    
                    
                    
            #if myListPrev[n]!=pileupread.qpos & pileupcolumn.n>10:
                    
            dataT=[]
            for i in range(len(readsA)):
                dataT.append(-1)
                    
            Matrix.append(myList)
            #print myList
            #if myList.count('A')+myList.count('C')+myList.count('T')+myList.count('G')>10 :
                #print "@COV",pileupcolumn.pos+1,myList.count('A'),myList.count('C'),myList.count('T'),myList.count('G'),myList.count('N')
            myList=[]
            for i in range(len(readsA)):
                myList.append('N')




#samfile.close()

data=array(data).transpose()

print data[0]
print "l=",len(data)
print "l=",len(data[0])





#print data[2035,:]
#print data[2036,:]
#print data[2037,:]


# computing K-Means with K = 2 (2 clusters)
centroids,_ = kmeans(data,2)
# assign each sample to a cluster
idx,_ = vq(data,centroids)





print idx
print "len(idx)",len(idx)



#list for allelese from cluster 0
a1=[]
a2=[]


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
            itemindex = np.where(readsA==pileupread.alignment.qname)
            #print "item=",itemindex[0]
            n=itemindex[0]
            #print n
            #print "PREVIOUS POSITION=", myListPrev[n]
            #print "current position=", pileupread.qpos
            
            if myListPrev[n]!=pileupread.qpos: # to skip splicing and deletions ; tehre is still base but it comes from previous
               
                #print idx[n]
                if pileupread.alignment.seq[pileupread.qpos].lower() != 'n' and idx[n]==0 :
                    a1.append(pileupread.alignment.seq[pileupread.qpos].lower())
                elif pileupread.alignment.seq[pileupread.qpos].lower() != 'n' and idx[n]==1 :
                    a2.append(pileupread.alignment.seq[pileupread.qpos].lower())
                myList=[]
                for i in range(len(readsA)):
                    myList.append('N')

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

                    
def printme( a1_a,a1_c,a1_t,a1_g):
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
    
    if idx[i]==0:
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

sys.exit(1)


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
    print '%s: %7d' % (letter, count)
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
    print i,letter
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
    print '%s: %7d' % (letter, count)
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


if len(haplo1)!=len(haplo2):
    print "ERROR!"
    print "exit"
    sys.exit(1)

print "len(haplo1)",len(haplo1)
print "len(haplo2)",len(haplo1)

if len(MatrixC1)>len(MatrixC2):
    ratio=float(len(MatrixC2))/len(MatrixC1)
else:
    ratio=len(MatrixC1)/float(len(MatrixC2))

print "ratio",ratio





for i in range(0, len(haplo1)):
    if i+g1+1>=g1E and i+g2+1<=g2E:
        if haplo1[i]!=haplo2[i] and haplo1[i]!='N' and haplo2[i]!='N':
            print "@SNP",i,i+g1+1,haplo1[i],haplo2[i], MatrixC1[:,i].tolist().count(haplo1[i]), MatrixC2[:,i].tolist().count(haplo2[i]),ratio
        elif haplo1[i]==haplo2[i] and haplo1[i]!='N' and haplo2[i]!='N':
            print "@HOMO",i,i+g1+1,haplo1[i],haplo2[i], MatrixC1[:,i].tolist().count(haplo1[i]), MatrixC2[:,i].tolist().count(haplo2[i]),ratio


print "MatrixC1 ", len(MatrixC1)
print "MatrixC2 ", len(MatrixC2)

print "NEW begin",g1
print "NEW end ",g2
print "ENSEMBL begin",g1E
print "ENSEMBL end ",g2E
print "ALLELE1", len(MatrixC1)
print "ALLELE2", len(MatrixC2)


haplo1E=[]
haplo2E=[]



#write haplotype staking into consideration coverage
for i in range(0, len(haplo1)):
    if i+g1+1>=g1E and i+g1+1<=g2E:
        print i+g1+1,g1E,g2E
        print MatrixC1[:,i].tolist().count(haplo1[i]), MatrixC2[:,i].tolist().count(haplo2[i])
        print haplo1[i]
        print haplo2[i]
        if  MatrixC1[:,i].tolist().count(haplo1[i])==1 and MatrixC2[:,i].tolist().count(haplo2[i])==1 :
            print "both N"
            haplo1E.append('N')
            haplo2E.append('N')
        elif MatrixC1[:,i].tolist().count(haplo1[i])==1:
            print "case1"
            haplo1E.append(haplo2[i])
            haplo2E.append(haplo2[i])
        elif MatrixC2[:,i].tolist().count(haplo2[i])==1:
            print "case2"
            haplo2E.append(haplo1[i])
            haplo1E.append(haplo1[i])
        else:
            print "Both good"
            haplo1E.append(haplo1[i])
            haplo2E.append(haplo2[i])
        if len(haplo1E)!=len(haplo2E):
            print haplo1E
            print haplo2E
            print "ERROR"
            sys.exit(1)






print "@HAPLO1",''.join(haplo1E)
print "@HAPLO2",''.join(haplo2E)

print "len(H1)",len(haplo1)
print "len(H2)",len(haplo2)
print "len(H1)",len(haplo1E)
print "len(H2)",len(haplo2E)


print "done!"


