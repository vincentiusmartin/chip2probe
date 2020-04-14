################## Parameters ###############
# Files and folders
protein = "runx1"
cell_line = 'jurkat'
genome_ver = 'hg38'
# -----------peakFile for tlymph------------------
peakFile = 'preprocessing/' + cell_line + '/' + protein + '_' + cell_line + '_peaks.bed'
# -----------peakFile for chipseq-----------------
#peakFile = 'preprocessing/chipseq/' + protein + '_idr_without_peaklist'
# -----------peakFile for dnaseq------------------
#peakFile = 'dnaseq_k562_peaks.narrowPeak'

#----------------genome file for hg19------------------
if genome_ver == 'hg19':
    genomeFile = 'preprocessing/genome/human_g1k_v37.fasta'
elif genome_ver == 'hg38':
    genomeFile = 'preprocessing/genome/hg38.fa'

kmerFile = 'preprocessing/' + protein + '_kmer_alignment.txt'
# kPosition of the core and where to center the call
# core is right exclusionary, left inclusive [)
if protein == 'ets1':
    core = (11,15) 
if protein == 'runx1':
    core = (12,17)
threshold = 0.38
isPalindrome = False
# Optional settings
saveDiagnostic = True # saves table with calling information for each seq
logFile = True
if protein == 'ets1':
    centerPos = 12
if protein == 'runx1':
    centerPos = 14
#############################################
if centerPos == 'default':
    centerPos = int((core[0] + core[1])/2)
##### Imports ####

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
import os
from pybedtools import BedTool
from Bio.Seq import reverse_complement, Seq

#### Functions ####
def convertToBed(df, isPalindrome):
    if isPalindrome == False:
        chrom, start, orient,scores = [],[],[],[]
        for row in zip(df['Chromosome'],df['Start'],df['centerPlus'],df['scorePlus'],df['centerMinus'],df['scoreMinus']):
            if row[2]:
                for centerP, score in zip(row[2], row[3]): # + sites
                    chrom.append(row[0])
                    start.append(row[1] + centerP)
                    orient.append('+')
                    scores.append(score)
            if row[4]:
                for centerN, score in zip(row[4],row[5]): # - sites
                    chrom.append(row[0])
                    start.append(row[1] + centerN)
                    orient.append('-')
                    scores.append(score)
    else:
        for row in zip(df['Chromosome'],df['Start'],df['centerPlus'],df['scorePlus']):
            if row[2]:
                for centerP, score in zip(row[2], row[3]): # + sites
                    chrom.append(row[0])
                    start.append(row[1] + centerP)
                    orient.append('+')
                    scores.append(score)
    bedDF = pd.DataFrame({'chrom':chrom,'start':start, 'end':start,'orient':orient, 'score':scores})
    bedDF['end'] = bedDF['end'] + 1 # exclusive end position
    bedDF = bedDF.drop_duplicates().sort_values(by=['chrom','start']) # some sites overlap, will call the same centers
    return(bedDF)

##### Read in kmer data and process ####
kmer = pd.read_csv(kmerFile, sep = '\t')
k = len(kmer['kmer'][0])
coreLen = core[1] - core[0]
# Find the kPositions required, any would be sufficient to call
if k > coreLen: 
    searchEnd = core[1]
    checkK = 0
    ReqKpos = set() #
    while checkK != core[0]:
        checkK = searchEnd - k
        if checkK <= core[0]:
            ReqKpos.add(checkK)
            searchEnd = searchEnd + 1
# Or find the group of all kPositions that are needed, all or none
else:
    searchStart = core[0]
    checkK = 0
    ReqKpos = set()
    while searchStart + k <= core[1]:
        ReqKpos.add(searchStart)
        searchStart = searchStart + 1
# Determine flanks of ReqKPos for threshold score reporting
ScoredKpos = ReqKpos.copy()
if k >= coreLen:
    ScoredKpos.add(min(ReqKpos) - 1)
    ScoredKpos.add(max(ReqKpos) + 1)        

# Generate dictionary for quick retreaval of relevant kmers
thrKmers = kmer[(kmer['Escore'] > threshold) & (kmer['kposition'].isin(ScoredKpos))]
kDict = dict(zip(thrKmers['kmer'],zip(thrKmers['kposition'],thrKmers['Escore'])))


def kmerMatch(seq):
    """
    Returns matched positions in the sequence and their kpositions
    Input: sequence and kmer length (k)
    Output: consecutive positions, kpositions, and scores above threshold
    """
    # Get the kposition and kscore for the peak, save a numpy array
    kpos,kscore = [], []
    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        if window in kDict:
            kpos.append(kDict[window][0])
            kscore.append(kDict[window][1])
        else:
            kpos.append(0)
            kscore.append(-0.5)
    kpos = np.array(kpos)
    kscore = np.array(kscore)
    # Get consecutive positions, kpositions, and score via numpy operations
    if k >= coreLen:
        position = list(filter(lambda x: len(x) != 1,np.split(np.r_[:len(kpos)], np.where(np.diff(kpos) != 1)[0]+1)))
        kpos = list(filter(lambda x: len(x) != 1,np.split(kpos, np.where(np.diff(kpos) != 1)[0]+1)))
    elif k < coreLen:
        reqLen = len(ReqKpos)
        position = list(filter(lambda x: len(x) == reqLen,np.split(np.r_[:len(kpos)], np.where(np.diff(kpos) != 1)[0]+1)))
        kpos = list(filter(lambda x: len(x) == reqLen,np.split(kpos, np.where(np.diff(kpos) != 1)[0]+1)))
    kScore = []
    for pos in position:
        kScore.append(kscore[pos])
    return(zip(position, kpos, kScore))

def findCenter(zippedMatch, orient, seqLen):
    """
    Given a zip of match position, kposition, and kscore
    Returns the center sites and threshold kscore
    """
    centerSites = []
    siteScores = []
    for pos, kpos, kScore in zippedMatch:
        centerSite = (centerPos - kpos[0]) + pos[0]
        if orient == 'rc':
            centerSite = (seqLen - centerSite) -1
        centerSites.append(centerSite)
        if k >= coreLen:
            score = threshold
            for score1, score2 in zip(kScore, kScore[1:]):
                caniScore = sorted([score1, score2])[0]
                if caniScore > score:
                    score = caniScore
            siteScores.append(score)
        elif k < coreLen:
            siteScores.append(min(kScore))
    return(pd.Series([centerSites, siteScores]))
             
#------------------------Prediction starts here-----------------------------------
chromosome = ""
start = 0
primer = "GTCTTGATTCGCTTGACGCTGCTG"
fwd = str(Seq("CGCGGTGCACTCTGGGAAATGTGGTTTTCGCGGCGC").reverse_complement()) + primer
#fwd = "CGCGGTGCACTCTGGGAAATGTGGTTTTCGCGGCGC"
fwd = "AGCGCAGCCGTAGACTCCGCTCAGATCCCCGGGTCGGTCTTGATTCGCTTGACGCTGCTG"
seq = Seq(fwd)
seq_len = len(fwd)
rev_comp = str(seq.reverse_complement())
#fwd = "CGCCAGCGCGGTGTGGTGGATCGTGATTCCAGCCCG"
col_names = ["Chromosome", "Start", "fwd", "seq_len", "rev_comp"]
lst = [[chromosome, start, fwd, seq_len, rev_comp]]
df = pd.DataFrame.from_records(lst, columns = col_names)

# Run kmerMatch on the peaks
if isPalindrome == True:
    df[["centerPlus","scorePlus"]] = df.apply(lambda df: findCenter(kmerMatch(df['fwd']), 'fwd', len(df['fwd'])), axis = 1) 
    df = df[~df['centerPlus'].isna()]
else:
    df[["centerPlus","scorePlus"]] = df.apply(lambda df: findCenter(kmerMatch(df['fwd']), 'fwd',len(df['fwd'])), axis = 1) 
    df[["centerMinus","scoreMinus"]] = df.apply(lambda df: findCenter(kmerMatch(df['rev_comp']), 'rc',len(df['rev_comp'])), axis = 1)
    df = df[~df['centerPlus'].isna() | ~df['centerMinus'].isna()]

finalBed = convertToBed(df, isPalindrome)

print(finalBed)

    
print("DOne!!")