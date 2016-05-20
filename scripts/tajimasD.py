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


# Genetic diversity parameters

def a1(n):
    '''Given n, this function returns the (n-1)th harmonic number'''
    n = int(n)
    return sum((1.0/d) for d in range(1,n))

def a2(n):
    '''Given n, this function returns the (n-1)th squared harmonic number'''
    n = int(n)
    return sum((1.0/(d**2)) for d in range(1,n))

def b1(n):
    '''Creates b1 for the variance of Tajima's theta'''
    n = int(n)
    return ((n+1) /  (3*(n-1)))

def b2(n):
    '''Creates b2 for the variance of Tajima's theta'''
    n = int(n)
    num = ((n**2) + n + 3) * 2
    den = 9 * n * (n-1)
    return num / den

def c1(n):
    n = int(n)
    '''Creates c1 for the variance of Tajima's theta'''
    return b1(n) - (1 / a1(n))

def c2(n):
    n = int(n)
    '''Creates c2 for the variance of Tajima's theta'''
    return b2(n) - ((n+2) / (a1(n) * n)) + (a2(n) / (a1(n) ** 2 ))

def e1(n):
    n = int(n)
    return c1(n) / a1(n)

def e2(n):
    n = int(n)
    return c2(n) / ( (a1(n) ** 2) + a2(n) )

def cn(n):
    n = int(n)
    num = ( 2 * n * a1(n) ) - (4 * (n-1))
    den = (n-1) * (n-2)
    return num / den

def vD(n):
    n = int(n)
    chunk1 = (a1(n) ** 2) / ( a2(n) + (a1(n) ** 2) )
    chunk2 = cn(n) - ((n + 1) / (n - 1))
    return 1 + (chunk1 * chunk2)

def uD(n):
    n = int(n)
    return a1(n) - 1 - vD(n)

# AKA # of segregating sites
def wattersons_theta(pop):
    '''Accepts a list containing tuples with the alignment'''
    K = len(pop)
    pop_size = len(pop[0])
    theta = K / a1(pop_size)
    return theta


def fu_and_li_theta(pop):
    count = 0
    for site in pop:
        if 1 in Counter(site).values():
            count += 1
    return count

def get_distance(seq_a, seq_b):
    diffs = 0
    length = len(seq_a)
    assert len(seq_a) == len(seq_b)
    for chr_a, chr_b in zip(seq_a, seq_b):
        if chr_a != chr_b:
            diffs += 1
    #return diffs / float(length)
    return diffs

def tajimas_theta(population):
    '''Accepts a list containing tuples wit the alignment'''
    pop_size = len(population[0])
    length = len(population)
    for i in range(length):
        haplotype_a = haplotypes[i]
        frequency_a = population[haplotype_a] / float(pop_size)
        for j in range(0, i):
            #for j in range(haplotype_count):
            haplotype_b = haplotypes[j]
            #frequency_b = population[haplotype_b] / float(pop_size)
            #frequency_pair = frequency_a * frequency_b
            #diversity += frequency_pair * get_distance(haplotype_a, haplotype_b)
            diversity += (get_distance(haplotype_a, haplotype_b))
    return (diversity * (2 / (pop_size * (pop_size-1))))

def fu_and_li_D(population):
    pop_size = len(population[0])
    S = len(population)
    a1_pop = a1(pop_size)
    num = S - ( a1_pop * fu_and_li_theta(population) )
    den = math.sqrt((uD(pop_size) * S) + (vD(pop_size) * (S ** 2)))
    return num / den


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
            theta_W = wattersons_theta(sample)
            theta_FL = fu_and_li_theta(sample)
            D_FL.append(fu_and_li_D(sample))
            # get Wattersons theta, pi, etc,
    return D_FL


def CV_KDE(oneD_array, expand = 1000):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=50) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    # add nothing to the end of grid and pdf so you can get a nice looking kde
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple

D_FL = np.asarray(aln2params())

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.scatter(DNAnumpy, CTCnumpy)
#ax.set_xlim([-1000,4000])
return_D_FL = CV_KDE(D_FL)

plt.plot(return_D_FL[0], return_D_FL[1],color = 'b', linestyle = '-', label="N = 1000, B = 1")
plt.xlabel('Fu and Lis D', fontsize = 18)
plt.ylabel('Probability', fontsize = 18)
output =  mydir + '/figs/test.png'
#plt.savefig(output)
plt.savefig(output, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.xscale()
plt.close()
