from __future__ import division
import os, argparse, random, math
from collections import Counter
from itertools import product
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
import pylab as plt

mydir = os.path.expanduser("~/github/GenomeDorm")


class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        '.fa' in file_lower:
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print "Not in FASTA format."

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list




def get_distance(seq_a, seq_b):
    diffs = 0
    length = len(seq_a)
    assert len(seq_a) == len(seq_b)
    for chr_a, chr_b in zip(seq_a, seq_b):
        if chr_a != chr_b:
            diffs += 1
    #return diffs / float(length)
    return diffs



def readAln(i):
    '''
    Takes in all the file path for the aligned fasta seqs, reads them using the classFASTA class,
    and outputs the alignment with a list of tuples, each tuple containing the
    alignment of a given site.
    '''
    class_test = classFASTA(i)
    ######### This command returns the nested list containing sequence names
    ######### and sequences to a flat list containing only sequences
    sequences = [ x[1] for x in class_test.readFASTA() ]
    zippedSeq = zip(*sequences)
    return zippedSeq


def aln2params():
    '''
    Reads all the alignments in a directory, returns a tuple containing the following:
    (N, Watterson's theta, pi, singletons)
    NOTE: Only looks at "silent"-site variation
    '''
    D_FL = []
    path = mydir + '/output_dir/pan_genome_sequences'
    for i in os.listdir(path):
        if i.endswith(".fa.aln"):
            name = i.split('.')[0]
            i = path + '/' + i
            align = readAln(i)
            list_i = []
            #print len(i)
            mod = len(align) %3
            if mod != 0:
                align = align[:-mod]
            # get every third
            every_third = align[2::3]
            S = 0

            sample = [x for x in every_third if (len(set(x)) > 1) and ('-' not in x)\
            and (len(x) ==11) ]
            if len(sample) ==0:
                continue
            print sample
            #theta_W = wattersons_theta(sample)
            #theta_FL = fu_and_li_theta(sample)
            #D_FL.append(fu_and_li_D(sample))
            # get Wattersons theta, pi, etc,
    return D_FL


def CV_KDE(oneD_array, expand = 1000):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    # add nothing to the end of grid and pdf so you can get a nice looking kde
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple

D_FL = np.asarray(aln2params())

#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#ax.scatter(DNAnumpy, CTCnumpy)
#ax.set_xlim([-1000,4000])
#return_D_FL = CV_KDE(D_FL)

#plt.plot(return_D_FL[0], return_D_FL[1],color = 'b', linestyle = '-', label="N = 1000, B = 1")
#plt.xlabel('Fu and Lis D', fontsize = 18)
#plt.ylabel('Probability', fontsize = 18)
#output =  mydir + '/figs/test.png'
#plt.savefig(output)
#plt.savefig(output, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.xscale()
#plt.close()
