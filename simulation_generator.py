# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 00:31:25 2015

@author: Harry Yang
"""

import numpy as np

def haplotype_simulator(num_pos,num_reads,ratio,num_alt_spl, diff_exp_ratio):
    """
    This function generates binary matrix of RNA Seq data
    IN : # POS for binary matrix(# cols), # reads (# rows), ratio(1:x where x = int)
    OUT : Binary matrix [M x N]
    """
    binary_matrix = np.zeros((num_reads,num_pos))
#    print binary_matrix
    number=np.random.random_integers(1,100)
#    print number
    haplo_basic=np.zeros((1,num_pos))
    haplo_one=np.zeros((1,num_pos))
    haplo_two=np.zeros((1,num_pos))
#   Make sample haplotypes
    for i in range(num_pos):
        number = np.random.random_integers(-1,0)
        if number == 0:
            haplo_basic[0][i] = 1
        elif number == -1:
            haplo_basic[0][i] = -1
        else:
            exit(1)
#    print haplo_basic
#   change haplo basic to two haplotypes
    for i in range(num_pos):
        number = np.random.randint(1,1000)
        if number % 20 == 1: # 20 = Prob of change - 5%
#            print i
            haplo_one[0][i]=1
            haplo_two[0][i]=-1
        else:
            haplo_one[0][i]=haplo_basic[0][i]
            haplo_two[0][i]=haplo_basic[0][i]
#   Test the percentage of change
#    var=0
#    for i in range(num_pos):
#        if haplo_one[0][i] != haplo_two[0][i]:
#            var=var+1
#    print var

#   exon splicing

    haplo_one_alt, haplo_two_alt = exon_splicing(haplo_one,haplo_two,num_alt_spl)
#    print haplo_one_alt, haplo_two_alt


#      MAKE BINARY MATRIX - DIFFERENT EXPRESSION PER HAPLO
#      Make an array that stores which row is from which haplotype
    origin_tracker = []
    for i in range(num_reads):
#          make two rand int - for which haplotype and which isoform
        haplo_marker = np.random.randint(1,1000)%diff_exp_ratio
        isoform_marker = np.random.randint(1,1000)%num_alt_spl
#          The diff_exp_ratio is in form of 1:diff_exp_ratio
#          if 0 -> then use haplo 2
        if haplo_marker == 0:
            binary_matrix[i]=haplo_two_alt[isoform_marker]
            origin_tracker.append(1)
        else:
            binary_matrix[i]=haplo_one_alt[isoform_marker]
            origin_tracker.append(0)
    return binary_matrix, origin_tracker




def exon_splicing(haplo_one,haplo_two,num_alt_spl):
    """
    IN : haplo_one(1 X n), haplo_two(1 X n), number of alternative splicing product
    OUT : haplo_one_alt(num_alt_spl X n), haplo_two_alt(num_alt_spl X n)
    """
    pos_size = len(haplo_one[0])
    haplo_one_alt=np.zeros((num_alt_spl,pos_size))
    haplo_two_alt=np.zeros((num_alt_spl,pos_size))
    spl_points=[]

    #SET SPLICING POINTS
    for i in range(num_alt_spl):
        rand_num=np.random.randint(1,pos_size-1)
        spl_points.append(rand_num)
    spl_points=sorted(spl_points)
#    print spl_points
    #COPY HAPLO TO ALT HAPLO
    haplo_one_alt[0]=haplo_one[0]
    haplo_two_alt[0]=haplo_two[0]

    #MAKE ALTERNATIVE SPLICING - Mutually exclusive exons
    for i in range(len(spl_points)-1):
        point_one=spl_points[i]
        point_two=spl_points[i+1]
        haplo_one_alt[i+1]=haplo_one_alt[0]
        haplo_two_alt[i+1]=haplo_two_alt[0]
        haplo_one_alt[i+1][point_one:point_two]=0
        haplo_two_alt[i+1][point_one:point_two]=0
    #POSSIBLY OTHER SPLICING STRUCTURES?
    #RETURN THE ALT SPLICINGS
    return haplo_one_alt, haplo_two_alt