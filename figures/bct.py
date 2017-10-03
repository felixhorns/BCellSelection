#!/usr/bin/env python

""" Functions for working with site frequency spectra (SFS) """

from __future__ import division
import sys
import os
import time
import copy
import itertools
import random

import numpy as np
import pandas as pd
import scipy

sys.path.append("/local10G/rfhorns/resources/betatree/src")
import betatree as bt
import sfs

from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns

# import logit scale
from LogitScale import *
mscale.register_scale(LogitScale)

def flatten(L):
    """ Convert list of lists into flat list"""
    return [item for sublist in L for item in sublist]

def unpack(X):
    """ Unpack a comma separated list of values into a flat list """
    return flatten([x.split(",") for x in list(X)])

def dict_to_tuple(D):
    """ Convert a dictionary using tuples as keys (K1, K2): V to a list of tuples (K1, K2, V)"""
    result = []
    for k, v in D.items():
        result.append(tuple(list(k) + [v]))
    return result

def get_muts_freqs(df_raw, lineage_uid):
    """ Get unique germline mutations and calculate their frequencies for a lineage """
    df = df_raw[df_raw["lineage_uid"] == lineage_uid]
    N = df.shape[0] # number of sequences in lineage
    df = df[df["mut_germline_positions"] != ""] # filter out sequences without germline mutations
    df = df.loc[-pd.isnull(df["mut_germline_positions"])] # filter out sequences without germline mutations
    muts = zip(unpack(df.mut_germline_positions), unpack(df.mut_germline_before), unpack(df.mut_germline_after))
    counts = {}
    for mut in muts:
        try:
            counts[mut] += 1
        except:
            counts[mut] = 1
    
    freqs = {}
    for k, v in counts.items():
        freqs[k] = float(v) / float(N)
    
    return freqs

def filter_VDJ_genes(freqs, keep=["V", "J"]):
    """ Removes mutations unless gene is in keep (e.g. V, J)"""
    freqs_filtered = {}
    for k, v in freqs.items():
        position, before, after = k
        if position[0] in keep:
            freqs_filtered[k] = v
    return freqs_filtered

def get_lineage_size(df_raw, lineage_uid):
    """ Get number of sequences in a lineage """
    df = df_raw[df_raw["lineage_uid"] == lineage_uid]
    N = df.shape[0] # number of sequences in lineage
    return N

def get_lineage_size_all(df):
    """ Get number of sequences for all lineages"""
    sizes = {}
    for lineage_uid in df.lineage_uid.unique():
        N = get_lineage_size(df, lineage_uid)
        sizes[lineage_uid] = N
    return sizes

def get_muts_freqs_all(df):
    """ Get germline mutations and calculate frequencies for all lineages"""
    
    start_time = time.time()
    
    freqs = {}
    for lineage_uid in df.lineage_uid.unique():
        myFreqs = get_muts_freqs(df, lineage_uid)
        myFreqs = filter_VDJ_genes(myFreqs)
        freqs[lineage_uid] = myFreqs

    elapsed_time = time.time() - start_time
    print "Wall clock time", elapsed_time
    print
    
    return freqs

def get_muts_counts(df_raw, lineage_uid):
    """ Get unique germline mutations and tally counts for a lineage """

    df = df_raw[df_raw["lineage_uid"] == lineage_uid]
    N = df.shape[0] # number of sequences in lineages
    
    df = df[df["mut_germline_positions"] != ""] # filter out sequences without germline mutations
    df = df.loc[-pd.isnull(df["mut_germline_positions"])] # filter out sequences without germline mutations

    muts = zip(unpack(df.mut_germline_positions), unpack(df.mut_germline_before), unpack(df.mut_germline_after))
    
    counts = {}
    for mut in muts:
        try:
            counts[mut] += 1
        except:
            counts[mut] = 1
            
    return counts, N

def get_muts_counts_all(df):
    """ Get germline mutations and tally counts for all lineages"""
    
    start_time = time.time()

    freqs = {}
    sizes = {}
    for lineage_uid in df.lineage_uid.unique():
        myFreqs, N = get_muts_counts(df, lineage_uid)
        myFreqs = filter_VDJ_genes(myFreqs)
        freqs[lineage_uid] = myFreqs
        sizes[lineage_uid] = N
        
    elapsed_time = time.time() - start_time
    print "Wall clock time", elapsed_time
    print
    
    return freqs, sizes

def load_germline_muts(indir):
    """ Load germline mutations for all subjects in a directory into a dictionary G{subject}{allele} = [mut, ...]"""
    
    G = {}
    
    for f in os.listdir(indir):
        if f.startswith("germline_muts_V6_Full_subject") and f.endswith(".txt"):
            
            patient_uid = f.split("germline_muts_V6_Full_subject")[1].replace(".txt", "")
            
            G[patient_uid] = {}
            
            with open(indir+f) as infile:
                for line in infile:
                    allele, mut = line.rstrip().split("\t")
                    try:
                        G[patient_uid][allele].append(mut)
                    except:
                        G[patient_uid][allele] = [mut]
                                    
    return G

def make_lineage_to_alleles(df):
    """ Get names of alleles that belong to each lineage """
    M = {}
    for lineage_uid, group in df.groupby("lineage_uid"):
        M[lineage_uid] = list(group["V_germline"].unique())
    return M

def drop_germline_muts(freqs, germline_muts, lineage_to_alleles):
    """ Drop germline mutations from mutation frequencies """
    
    verbose = False
    
    freqs_clean = {}
    
    dropped_count = 0
    
    for lineage_uid, my_muts in freqs.items():
        
        # Get patient uid and alleles of lineage
        my_patient_uid = str(lineage_uid)[0]
        my_alleles = lineage_to_alleles[lineage_uid]
        
        # Look up germline mutations for this patient and allele
        my_germline_muts = []
        for x in my_alleles:
            if x in germline_muts[my_patient_uid].keys():
                my_germline_muts.extend(germline_muts[my_patient_uid][x])
        my_germline_muts = list(set(my_germline_muts))
        
        # Walk through mutations and keep mutations that are not germline mutations
        for key, freq in my_muts.items():
            
            position, before, after = key
            mut_str = position + "_" + before + "_" + after
            
            if mut_str in my_germline_muts and verbose:
                print mut_str, freq
            
            if mut_str not in my_germline_muts:
                try:
                    freqs_clean[lineage_uid][key] = freq
                except:
                    freqs_clean[lineage_uid] = {key: freq}
            else:
                dropped_count += 1

    print "Dropped germline mutations", dropped_count
    return freqs_clean

def bin_sfs(freqs, bins):
    binned_sfs, _ = np.histogram(freqs, bins=bins) # count mutations per bin
    bin_sizes = bins[1:] - bins[:-1]
    binned_sfs_normed = binned_sfs / bin_sizes # normalize by size of bin
    return binned_sfs, binned_sfs_normed

def bin_sfs_cut(freqs, bins, leaves):
    # cuts SFS at appropriate size for the lineage by setting all bins beyond 1/leaves to np.nan
    binned_sfs, _ = np.histogram(freqs, bins=bins) # count mutations per bin
    bin_sizes = bins[1:] - bins[:-1]
    binned_sfs_normed = binned_sfs / bin_sizes # normalize by size of bin
    binned_sfs_normed[bins[1:] < 1/leaves] = np.nan
    binned_sfs_normed[bins[:-1] > 1-(1/leaves)] = np.nan
    return binned_sfs, binned_sfs_normed

def unpack_freqs(D):
    """ Unpacks a dictionary {lineage_uid: {mut: freq, ...}}
        to a list of freqs """
    freqs = flatten([x.values() for x in D.values()])
    return freqs

def drop_zeroes(X, Y):
    """ Drops data points where y = 0 """
    x_dropped = []
    y_dropped = []
    for x, y in zip(X, Y):
        if y != 0:
            x_dropped.append(x)
            y_dropped.append(y)
    return x_dropped, y_dropped

def plot_sfs(ax, bin_centers, binned_sfs, p_min=1e-4, **kwargs):
    X_clean, Y_clean = drop_zeroes(bin_centers, binned_sfs)
    ax.plot(X_clean, Y_clean, **kwargs)
    ax.set_xscale("logit", p_min=p_min)
    ax.set_yscale("log")
    ax.set_xlabel("Mutation frequency")
    ax.set_ylabel("Density of mutations")
    sns.despine()
    plt.tight_layout()
    return ax

def calc_average_sfs(freqs, lineage_sizes, lineage_uids, bins):
    """ Calculate mean SFS over an ensemble by taking mean value at each bin """
    S = np.empty((len(lineage_uids), len(bins)-1))
    for i, lineage_uid in enumerate(lineage_uids):
        myFreqs = freqs[lineage_uid].values()
        myLeaves = lineage_sizes[lineage_uid]
        binned_sfs, binned_sfs_normed = bin_sfs_cut(myFreqs, bins=bins, leaves=myLeaves)
        S[i,:] = binned_sfs_normed
    # S[S==0.0] = np.nan # before calculating mean, convert 0 to nan
    # S_mean = np.nanmedian(S, axis=0)
    S_mean = np.nanmean(S, axis=0)    
    return S_mean

# Simulate ensemble of SFSs
def simSFSs(alpha=2, n=1000, ntrees=1, N=1, mode = "log", bins=None):
    if bins is None: bins = np.array([1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999])
    start_time = time.time()
    SFSs = []
    for _ in xrange(N):
        mySFS = sfs.SFS(n,alpha=alpha)
        mySFS.getSFS(ntrees = ntrees)
        mySFS.binSFS(mode=mode, bins=bins) # bin SFS
        SFSs.append(mySFS)
    print "Elapsed time (wall clock):", time.time() - start_time
    return SFSs

# def simSFSs_sampleLineageSizes(pop_sizes, alpha=2, ntrees=1, bins=bins):
#     """ Simulate evolution with population sizes sampled from distribution and calculate SFS """
#     start_time = time.time()
#     SFSs = []
#     for pop_size in pop_sizes:
#         mySFS = sfs.SFS(pop_size,alpha=alpha)
#         mySFS.getSFS(ntrees=ntrees)
#         mySFS.binSFS(bins=bins) # bin SFS
#         SFSs.append(mySFS)
#     print "Elapsed time (wall clock):", time.time() - start_time
#     return SFSs

def simSFSs_sampleLineageSizes(pop_sizes, alpha=2, ntrees=1, bins=None):
    """ Simulate evolution with population sizes sampled from distribution and calculate SFS """
    if bins is None: bins = np.array([1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999])
    start_time = time.time()
    SFSs = []
    for _ in xrange(len(pop_sizes)):
        pop_size = random.choice(pop_sizes)
        mySFS = sfs.SFS(pop_size,alpha=alpha)
        mySFS.getSFS(ntrees=ntrees)
        mySFS.binSFS(bins=bins) # bin SFS
        SFSs.append(mySFS)
    print "Elapsed time (wall clock):", time.time() - start_time
    return SFSs

def censor_binned_sfs(binned_sfs, leaves, bins=None):
    """ Censor bins that cannot have a value based on population size """
    if bins is None: bins = np.array([1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999])
    sfs_censored = copy.deepcopy(binned_sfs)
    sfs_censored[bins[1:] <= 1/leaves] = np.nan
    sfs_censored[bins[:-1] >= 1-(1/leaves)] = np.nan
#     sfs_censored[bins[1:] <= 1/leaves] = 0.
#     sfs_censored[bins[:-1] >= 1-(1/leaves)] = 0.
    return sfs_censored

def censor_mean_SFSs(SFSs, lineage_sizes, bins=None):
    """ Censor SFSs, compute mean and SEM """
    if bins is None: bins = np.array([1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999])
    S = np.zeros((len(lineage_sizes), len(bins)-1))
    for i, (SFS, lineage_size) in enumerate(zip(SFSs, lineage_sizes)):
        S[i,:] = censor_binned_sfs(SFS.binned_sfs, lineage_size, bins)
    S_mean = np.nanmean(S, axis=0)
#    S_mean = np.nanmedian(S, axis=0)
    S_sem = scipy.stats.sem(S, axis=0)
    return S_mean, S_sem

def drop_zeroes_errors(X, Y, Z):
    """ Drops data points where y = 0 or z = 0"""
    x_dropped = []
    y_dropped = []
    z_dropped = []
    for x, y, z in zip(X, Y, Z):
        if y != 0 or z != 0:
            x_dropped.append(x)
            y_dropped.append(y)
            z_dropped.append(z)
    return x_dropped, y_dropped, z_dropped

