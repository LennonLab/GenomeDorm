from __future__ import division
import os, argparse, random, math
from collections import Counter
from itertools import product, groupby
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
import pylab as plt
import pandas as pd

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

class popGenStats:

    def __init__(self, pop):
        self.pop = pop
        self.n = len(self.pop[0])
        self.S = len(self.pop)

    def a1(self):
        '''Given n, this function returns the (n-1)th harmonic number'''
        return sum((1.0/d) for d in range(1,self.n))

    def a2(self):
        '''Given n, this function returns the (n-1)th squared harmonic number'''
        #n = int(n)
        return sum((1.0/(d**2)) for d in range(1,self.n))

    def b1(self):
        '''Creates b1 for the variance of Tajima's theta'''
        #n = int(n)
        return ((self.n+1) /  (3*(self.n-1)))

    def b2(self):
        '''Creates b2 for the variance of Tajima's theta'''
        #n = int(n)
        num = ((self.n**2) + self.n + 3) * 2
        den = 9 * self.n * (self.n-1)
        return num / den

    def c1(self):
        #n = int(n)
        '''Creates c1 for the variance of Tajima's theta'''
        return self.b1() - (1 / self.a1())

    def c2(self):
        #n = int(n)
        '''Creates c2 for the variance of Tajima's theta'''
        return self.b2() - ((self.n+2) / (self.a1() * self.n)) + (self.a2() / (self.a1() ** 2 ))

    def e1(self):
        #n = int(n)
        return self.c1() / self.a1()

    def e2(self):
        #n = int(n)
        return self.c2() / ( (self.a1() ** 2) + self.a2() )

    def cn(self):
        #n = int(n)
        num = ( 2 * self.n * self.a1() ) - (4 * (self.n-1))
        den = (self.n-1) * (self.n-2)
        return num / den

    def vD(self):
        #n = int(n)
        chunk1 = (self.a1() ** 2) / ( self.a2() + (self.a1() ** 2) )
        chunk2 = self.cn() - ((self.n + 1) / (self.n - 1))
        return 1 + (chunk1 * chunk2)

    def uD(self):
        #n = int(n)
        return self.a1() - 1 - self.vD()

    # AKA # of segregating sites
    def wattersons_theta(self):
        '''Accepts a list containing tuples with the alignment'''
        theta = self.S / self.a1()
        return theta


    def fu_and_li_theta(self):
        count = 0
        for site in self.pop:
            if 1 in Counter(site).values():
                count += 1
        return count


    def tajimas_theta(self):
        '''Accepts a list containing tuples wit the alignment'''
        #n = len(population[0])
        pi_sum = 0
        for site in self.pop:
            #n = len(site)
            site_dict = Counter(site)
            sum_freq = 0
            for key, value in site_dict.items():
                freq = value / self.n
                site_dict[key] = freq
                sum_freq += freq ** 2
            pi_sum += (self.n/(self.n-1.0)) * (1 - sum_freq)
        return pi_sum

    def tajimas_D(self):
        num = self.tajimas_theta() - self.wattersons_theta()
        den = math.sqrt( (self.e1() * self.S) + (self.e2() * self.S * (self.S-1)) )
        T_D = num / den
        return T_D

    def fu_and_li_D(self):
        num = self.S - ( self.a1() * self.fu_and_li_theta() )
        den = math.sqrt((self.uD() * self.S) + (self.vD() * (self.S ** 2)))
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
    align_dict = {}
    path = mydir + '/data/out_test_1465341694/pan_genome_sequences'
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
            and (len(x) ==14) ]
            if len(x) != 14 and len(x) > 1:
                print name, str(len(x))

            if len(sample) ==0:
                continue
            testtt =  popGenStats(sample)
            pi = popGenStats(sample).tajimas_theta()
            theta_W = popGenStats(sample).wattersons_theta()
            theta_FL = popGenStats(sample).fu_and_li_theta()
            T_D = popGenStats(sample).tajimas_D()
            D_FL = popGenStats(sample).fu_and_li_D()
            align_dict[name] = [pi, theta_W, theta_FL, T_D, D_FL]
            # get Wattersons theta, pi, etc,
    df = pd.DataFrame(align_dict, index =None)
    df = df.transpose()
    #df.index = range( df.shape[0])
    df = df.reset_index()
    #print df.shape[0]
    print df
    df.columns = ['Gene', 'W_T', 'pi', 'FL_T', 'T_D', 'FL_D']
    df.to_csv(mydir + '/data/out/PopGenStats.txt', sep = '\t', index=False)
    return df


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



def selectDormGenes():
    dataframe = pd.read_csv(mydir + '/data/out/PopGenStats.txt', sep = '\t')
    PresAbs = pd.read_csv(mydir + '/data/out_test_1465341694/gene_presence_absence.csv')
    spore_genes = PresAbs[PresAbs['Annotation'].str.contains("spor")]
    spore_genes_names = spore_genes[['Gene']]
    print spore_genes_names
    #print spore_genes_names
    print pd.merge(dataframe, spore_genes_names, on=['Gene'])

#selectDormGenes()


def makePlot(param = 'pi'):
    dataframe = pd.read_csv(mydir + '/data/out/PopGenStats.txt', sep = '\t')
    param_np = np.asarray(dataframe[[param]])
    getKDE = CV_KDE(param_np)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.plot(getKDE[0], getKDE[1],color = 'b', linestyle = '-', label="N = 1000, B = 1")
    if param == 'T_D':
        plt.xlabel('Tajimas D', fontsize = 18)
    elif param == 'FL_D':
        plt.xlabel('Fu and Lis D', fontsize = 18)
    elif param == 'W_T':
        plt.xlabel('Wattersons Theta', fontsize = 18)
    else:
        plt.xlabel('pi', fontsize = 18)
    plt.ylabel('Probability', fontsize = 18)
    output =  mydir + '/figs/' + str(param) + '.png'
    plt.savefig(output, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

params = ['pi', 'T_D', 'FL_D','W_T']
for param in params:
     makePlot(param = param)
