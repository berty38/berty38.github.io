# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 13:58:24 2018

A tool set to help discover "counterfactual" subsequences in a time-series.

@author: colin shea-blymyer
"""

import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from multiprocessing import cpu_count
from tqdm import trange, tqdm


def sequence_scores(data, init_1, init_2, window_size):
    """ Return 2-list of indicator arrays for each initial condition.
    
        Arguments:
        data        -- the data set to analyze
        init_1      -- the first desired initial condition
        init_2      -- the second desired initial condition
        window_size -- the desired size of each subsequence
        
        Returns:
        A 2-item list. The first item is a list that contains an entry for each
        value in data that equals init_1, as long as that value is not within
        window_size values of the end of data. Each entry indicates the 
        negative of the number of values equal to either init_1 or init_2 
        within window_size values to the entry's right. The second item does
        the same for init_2.
    """
    offset = window_size - 1
    scores_1 = np.zeros(len(data)-offset)
    scores_2 = np.zeros(len(data)-offset)
    for i, x in enumerate(data[0:-offset]):
        scan_size = window_size
        if i >= (len(data)-window_size):
            scan_size = (len(data)-(window_size-1)-i)+1
        if x == init_1:
            scores_2[i] = None
            for s in data[i+1: i+scan_size]:
                if s == init_1 or s == init_2:
                    scores_1[i] -= 1
        elif x == init_2:
            scores_1[i] = None
            for s in data[i+1: i+scan_size]:
                if s == init_1 or s == init_2:
                    scores_2[i] -= 1
        else:
            scores_1[i] = None
            scores_2[i] = None
                    
    return [scores_1, scores_2]
    

def start_select(score_1, score_2, window_size):
    """ Return a 2-tuple of subsequence starting indeces.
    
        Arguments:
        score_1     -- the first score as returned by sequence_scores
        score_2     -- the second score as returned by sequence_scores
        window_size -- the desired size of each subsequence
        
        Returns:
        A 2-tuple. The first item is a list of indeces that indicate the start
        of a desired subsequence that corresponds with score_1. The second item
        is as such for score_2.
        
        Notes:
        This method begins with the score list that has the largest "limit
        rate" - the ratio of the sum of the array to the number of scores in
        the array. It then chooses the index of the score with the highest
        score, making all other scores within window_size indeces right
        of it invalid for future selection. The method then does the same for
        the score list with the smaller limit rate. This process is repeated
        until there are no more valid scores available, checking at each step
        to ensure that the score selected does not contain previously selected
        scores.
    """
    # pad the ending of the score arrays
    for i in range(window_size-1):
        score_1 = np.append(score_1, np.nan)
        score_2 = np.append(score_2, np.nan)
        
    limit_rat_1 = limit_rate(score_1)
    limit_rat_2 = limit_rate(score_2)
    
    # lists of indices that start selected subsequences
    seq_1 = []
    seq_2 = []    
    if limit_rat_1 <= limit_rat_2:
        while np.nanmax(score_1) != np.NINF and np.nanmax(score_2) != np.NINF:
            score_1, score_2, seq_1 = update_selection(score_1, score_2, seq_1, window_size, -1)
            score_2, score_1, seq_2 = update_selection(score_2, score_1, seq_2, window_size, -1)                
                
    else:
        while np.nanmax(score_1) != np.NINF and np.nanmax(score_2) != np.NINF:
            score_2, score_1, seq_2 = update_selection(score_2, score_1, seq_2, window_size, -1)
            score_1, score_2, seq_1 = update_selection(score_1, score_2, seq_1, window_size, -1)
                
    return (seq_1, seq_2)
            

def update_selection(score_1, score_2, seq_1, window_size, index):
    """Helper funciton to check subsequences for availability."""    
    
    # find index for subsequence starting from last high score
    i = int(np.where(score_1 == np.nanmax(score_1))[-1][index])
    
    # check if sequence is in use
    if not np.NINF in score_1[i:i+window_size]:
        # add index to sequences, and show subsequence in use for both            
        seq_1.append(i)
        score_1[i:i+window_size] = np.NINF
        score_2[i:i+window_size] = np.NINF
    else:
        # sequence is in use, this starting point, and any others
        # contained within the possible subsequence are invalid,
        new_sub_1 = score_1[i:i+window_size]
        new_sub_2 = score_2[i:i+window_size]
        new_sub_1[np.where(new_sub_1 != np.NINF)] = np.nan
        new_sub_2[np.where(new_sub_2 != np.NINF)] = np.nan
        score_1[i:i+window_size] = new_sub_1
        score_2[i:i+window_size] = new_sub_2
        # try again if possible
        if np.nanmax(score_1) != np.NINF:
            score_1, score_2, seq_1 = update_selection(score_1, score_2, seq_1, window_size, index)
    
    return score_1, score_2, seq_1


def limit_rate(score):
    """Compute the ratio the array's sum to its number of scores."""
    num_valid = np.count_nonzero(~np.isnan(score))
    score_sum = np.nansum(score)
    limit_rate = score_sum/num_valid
    return limit_rate


def select_sequences(indeces, data, window_size):
    """Returns the subsequences of window_size length starting at indeces"""
    subsequences = []    
    for index in indeces:
        subsequence = data[index:index+window_size]
        subsequences.append(subsequence)
    
    return subsequences


def score_and_start(data, init1, init2, window_size):
    """ Returns a 2-tuple of indeces.
    
        The entries of the 2-tuple correspond to init1 and init2, respectively.
        Each value in each entry corresponds to the start of a subsequence of 
        data with lenth equal to window_size.        
    """
    scores = sequence_scores(data, init1, init2, window_size)
    return start_select(scores[0], scores[1], window_size)


def monte_carlo_profile(depth=5):
    """ Profiles how this method compares to random subsequence selection."""
    # how does algorithm perform with random data, with different parameters?
    # compared to random selection?
    # create distance matrix of difference in performance of random and select
    num_cores = cpu_count() - 2
    
    diffs = np.zeros((depth, depth))
    params = []
    data_sizes = np.linspace(1000, 1000000, num=depth).astype(int)
    for data_size in data_sizes:
        window_sizes = np.linspace(np.log(data_size), np.sqrt(data_size), num=depth).astype(int)
        for window_size in window_sizes:
            data = np.random.randint(0, 100, data_size)
            params.append((data, window_size))
        
    random_scores_list = []
    for i in trange(5):
        random_selections = Parallel(n_jobs=num_cores)(delayed(select_random)(param[0], 25, 75, param[1]) for param in params)
        random_scores = Parallel(n_jobs=num_cores)(delayed(performance_score)(selection) for selection in random_selections)
        random_scores_list.append(random_scores)
        
    random_scores_list = np.asarray(random_scores_list)
    random_scores_avg = np.mean(random_scores_list, axis=0)
    
    selections = Parallel(n_jobs=num_cores)(delayed(score_and_start)(param[0], 25, 75, param[1]) for param in tqdm(params))
    select_scores = Parallel(n_jobs=num_cores)(delayed(performance_score)(selection) for selection in selections)        
    
    diffs = select_scores - random_scores_avg
    diffs = diffs.reshape((depth, depth))
            
    fig, ax = plt.subplots()
    im = ax.imshow(diffs)
    ax.set_yticks(np.arange(len(data_sizes)))
    ax.set_yticklabels(data_sizes)
    
    ax.set_xticks(np.arange(len(data_sizes)))
    
    for i in range(depth):
        for j in range(depth):
            text = ax.text(j, i, diffs[i, j], ha="center", va="center", color="w")
    
    ax.set_title("Difference between number of subsequences discovered by heuristic and random selection")
    plt.xlabel("Window Size - linear spacing from log of data size to square root of data size")
    plt.ylabel("Data Size")
    fig.tight_layout()
    plt.show()

    

def performance_score(results):
    """ Scores the method, given by the number of sequences it finds."""
    scores = []
    scores.append(len(results[0]))
    scores.append(len(results[1]))
    return min(scores)
    

def select_random(data, init1, init2, window_size):
    """ Return a 2-tuple of subsequence starting indeces.
    
        Arguments:
        data        -- the data set to analyze
        init_1      -- the first desired initial condition
        init_2      -- the second desired initial condition
        window_size -- the desired size of each subsequence
        
        Returns:
        A 2-tuple. The first item is a list of indeces that indicate the start
        of a desired subsequence that corresponds with init1. The second item
        is as such for init2.
        
        Notes:
        This method selects random subsequences until no more subsequences can
        be selected.
    """
    seq_1 = []
    seq_2 = []
    i1 = np.where(data[0:-(window_size-1)] == init1)[0]
    i2 = np.where(data[0:-(window_size-1)] == init2)[0]
    while i1.any() and i2.any():
        s1 = np.random.choice(i1)
        seq_1.append(s1)
        i1 = np.delete(i1, np.where((i1>=s1-window_size) & (i1<=s1+window_size)))
        i2 = np.delete(i2, np.where((i2>=s1-window_size) & (i2<=s1+window_size)))
        if i2.any():
            s2 = np.random.choice(i2)
            seq_2.append(s2)
            i2 = np.delete(i2, np.where((i2>=s2-window_size) & (i2<=s2+window_size)))
            i1 = np.delete(i1, np.where((i1>=s2-window_size) & (i1<=s2+window_size)))
            
    return (seq_1, seq_2)
