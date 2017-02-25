# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 00:11:35 2015

@author: Harry Yang
"""

import sklearn as learn
from sklearn import cluster
import numpy as np
import scipy
import pandas
from simulation_generator import haplotype_simulator
import os
import sys

def similarity(matrix):
    """
    compare similarity of two matrix.
    input: 2xlen matrix where len = length of haplo
    output: sim. score (percentage)
    """
    sim_score=0; zero_score=0
    non_zero_counter=0
    num_pos=len(matrix[0])
    for i in range(len(matrix[0])):

        if matrix[0][i] == 0:
            zero_score = zero_score+1
        elif matrix[1][i] == 0:
            zero_score = zero_score+1
        elif matrix[0][i]==matrix[1][i]:
            sim_score=sim_score+1
            non_zero_counter=non_zero_counter+1
        elif matrix[0][i] != matrix[1][i]:
            non_zero_counter=non_zero_counter+1
#    return sim_score*1.0/(zero_score*1.0/len(matrix[0]))*non_zero_counter, zero_score*1.0/len(matrix[0])
    if non_zero_counter == 0:
        non_zero_counter = 0
    return (sim_score*1.0/non_zero_counter)*(non_zero_counter), zero_score*1.0/len(matrix[0])


def test(repeat):
    sim_avg=0.0;
    num_successful_clustering=0;
    for i in range(repeat):
        test_matrix,origin_tracker=haplotype_simulator(1000,1000,100,2,3)
        print origin_tracker
        test_center=[]
        test_center.append(test_matrix[0])
        test_center.append(test_matrix[500])
        sim_avg = sim_avg + similarity(test_center)
#        test_bandwidth=learn.cluster.estimate_bandwidth(test_matrix, quantile=0.5)
        ms_vector = learn.cluster.MeanShift(seeds=test_center).fit_predict(test_matrix)
        if ms_vector[0] == ms_vector[251] and ms_vector[523]==ms_vector[809]:
            num_successful_clustering=num_successful_clustering+1
    print "average simliarty = ", sim_avg/repeat, "and # successful clustering = ", num_successful_clustering

def mean_shift_clustering(matrix):

    print len(matrix)
    center, sim_score, center_one_row, center_two_row = find_minimum_pair_new(matrix,len(matrix),0)
    print np.count_nonzero(center[0]), np.count_nonzero(center[1]), sim_score
    ms_vector = learn.cluster.MeanShift(seeds=center).fit_predict(matrix)
    return ms_vector

def test_complex(repeat, diff_ratio, num_trial, size):
    """
    diff_ratio used for simulator
    num_trial = num_trials to find center - use minimum pair
    """
    sim_avg=0.0
    num_successful_clustering=0
    cluster_score_avg=0.0
    inaccurate_min_pair_counter = 0
    for i in range(repeat):
        test_matrix,origin_tracker=haplotype_simulator(1000,size,100,3,diff_ratio)
#        find_minimum_pair_new(test_matrix, size, num_trial)
        test_center, sim_score, center_one_row, center_two_row = find_minimum_pair_new(test_matrix, size, num_trial)
        if origin_tracker[center_one_row] == origin_tracker[center_two_row]:
            inaccurate_min_pair_counter = inaccurate_min_pair_counter + 1
        sim_avg = sim_avg + sim_score
        print
        ms_vector = learn.cluster.MeanShift(seeds=test_center).fit_predict(test_matrix)
        cluster_score = cluster_checker(origin_tracker,ms_vector)
        if cluster_score == 0.0:
            num_successful_clustering = num_successful_clustering + 1

        pair_ind=""
        if origin_tracker[center_one_row] == origin_tracker[center_two_row]:
            pair_ind = "WRONG PAIR"
        else:
            pair_ind = "RIGHT PAIR"
        print "Trial # ", i, "Error Percentage ",cluster_score, "n(NONZERO) ", np.count_nonzero(test_center[0]), np.count_nonzero(test_center[1]), pair_ind
        cluster_score_avg = cluster_score_avg + cluster_score



    print "average simliarty = ", sim_avg/repeat, "and # successful clustering = ", num_successful_clustering
    print "average success rate = ", 1-(cluster_score_avg/repeat), "# wrong min pair = ", inaccurate_min_pair_counter


def cluster_checker(origin_tracker, cluster_vector):
    if len(origin_tracker) != len(cluster_vector):
        print "Error - the size of each vector is different"
    inaccuracy_tracker=0
    ori_haplo_one=0; ori_haplo_two=1
    cls_haplo_one=-1; cls_haplo_two=-1
#    check first one for comparison
    if origin_tracker[0] == 0 and cluster_vector[0] == 0:
        cls_haplo_one = 0; cls_haplo_two=1
    elif origin_tracker[0] == 0 and cluster_vector[0] == 1:
        cls_haplo_one = 1; cls_haplo_two=0
    elif origin_tracker[0] == 1 and cluster_vector[0] == 0:
        cls_haplo_one = 1; cls_haplo_two=0
    elif origin_tracker[0] == 1 and cluster_vector[0] == 1:
        cls_haplo_one = 0; cls_haplo_two=1

    for i in range(len(origin_tracker)):
#       origin_tracker is in either 1 or 0 and so is cluster vector
#       However, it is ambiguous whether 0 in origin_tracker corresponds to 0 in cluster_vector
        if origin_tracker[i] == ori_haplo_one and cluster_vector[i] != cls_haplo_one:
            inaccuracy_tracker=inaccuracy_tracker+1
        if origin_tracker[i] == ori_haplo_two and cluster_vector[i] != cls_haplo_two:
            inaccuracy_tracker=inaccuracy_tracker+1
    return inaccuracy_tracker*1.0/len(origin_tracker)



def find_minimum_pair(matrix, size, num_trial):
    """
    find minimum pair in the matrix randomly
    """
    center=[]
    center_one_row=0; center_two_row=0;
    sim_score=1
    zero_score_factor=0.9
    for i in range(num_trial):
        marker_one=np.random.randint(1,size)
        marker_two=np.random.randint(1,size)
        dummy_matrix=[]
        dummy_matrix.append(matrix[marker_one])
        dummy_matrix.append(matrix[marker_two])
        dummy_sim_score_pre, dummy_zero_score = similarity(dummy_matrix)
        dummy_sim_score = dummy_sim_score_pre + (zero_score_factor*dummy_zero_score)
        if dummy_sim_score < sim_score :
            center_one_row = marker_one
            center_two_row = marker_two
            sim_score = dummy_sim_score
    center.append(matrix[center_one_row])
    center.append(matrix[center_two_row])
    return center, sim_score, center_one_row, center_two_row

def find_minimum_pair_new(matrix,size,num_trial):
    center = []
    center_one_row=-1; center_two_row=-1;
    sim_score=np.inf
    zero_score_factor=1.0
#    print matrix.shape[0]
    dummy_matrix=[]
    nonzero_counter=0.90*np.count_nonzero(matrix[0])
#    print "NONZERO_COUNTER", nonzero_counter
    num_reads=0
    origin_matrix=[]
#    get the number of non_zero elements and add reads that are higher than the threshold
    for i in range(size-1):
#        print np.count_nonzero(matrix[i])
        if np.count_nonzero(matrix[i]) > nonzero_counter:
#            nonzero_counter=np.count_nonzero(matrix[i])
            dummy_matrix.append(matrix[i])
            num_reads=num_reads+1
            origin_matrix.append(i)
#            print "ADDED"
#    print "NUM_READS IN CRITERIA",num_reads


#    find best pair
    for j in range(num_reads):
        for k in range(j,num_reads):
            dummy_two_matrix=[]
            dummy_two_matrix.append(dummy_matrix[j])
            dummy_two_matrix.append(dummy_matrix[k])
            dummy_nonzero_num = np.count_nonzero(dummy_matrix[k])
            if 0.8*dummy_nonzero_num < np.count_nonzero(dummy_matrix[j]):
                if np.count_nonzero(dummy_matrix[j]) <= 1.5*dummy_nonzero_num:
                    dummy_sim_score, dummy_zero = similarity(dummy_two_matrix)
                    if dummy_sim_score < sim_score:
                        sim_score = dummy_sim_score
                        center_one_row = origin_matrix[j]
                        center_two_row = origin_matrix[k]
                else:
                    continue
                    # sys.exit(WRONG_NUMBER_READS)
            else:
                continue
                # sys.exit(WRONG_NUMBER_READS)
    if center_one_row == -1:
        sys.exit(25)
    # print "SEED MATRIX", matrix[center_one_row], matrix[center_two_row]
    center.append(matrix[center_one_row])
    center.append(matrix[center_two_row])
    return center, sim_score, center_one_row, center_two_row










if __name__ == "__main__":
#    print "hello world"
#    t= np.zeros((3,3), dtype = np.double)
#    t[0:1,0]=1
#    t[0:1,1]=0
#    t[0:2.2]=1
#    test=test_matrix_40_reads_100_postxt
#    for i in range(1,101):
#        test[0,i]=test[1,i]
#    k=cluster.k_means(test,n_clusters=2)
#    print k[1]
#
#    for i in range(len(k[1])):
#        print k[1][i]
#    test(10)

     test_complex(10,2,300,100)
#    print cluster_checker([1, 0, 0, 1, 0, 0], [1, 0, 1, 1,1, 0])
#    print cluster.k_means(test_matrix,n_init=1000,n_clusters=2,max_iter=1000,precompute_distances=True)[1]
#    print cluster.AgglomerativeClustering(n_clusters=2, affinity='euclidean').fit_predict(test_matrix)
#    print cluster.spectral_clustering(test_matrix, n_clusters=2)