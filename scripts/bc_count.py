#Load packages
import re
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os

#UPDATE ME!
reffile = 'unit_test_refs.fasta'
fastqdir = 'fastqs/A'
refstrain = 'ref'

rep = os.path.basename(fastqdir)
os.makedirs('bc_count_results', exist_ok=True)

#Build list of bc reference sequences
print('making list of sequences')
bcrefs = SeqIO.parse(reffile, "fasta")
bc_seqs = []
seqnames = []
for ref in bcrefs:
    bc_seqs.append(str(ref.seq))
    seqnames.append((str(ref.seq),ref.id))

namesdf = pd.DataFrame(seqnames, columns=['seq', 'strain'])

alldfs = []

#prep
print('counting')
for file in os.listdir(fastqdir):
    if not (file.endswith('.fastq')):
        print(file)
        continue
    filepath = os.path.join(fastqdir, file)
    basename, extension = os.path.splitext(file)
    
    allseqs = []
    bccounts = []
    
#list of sequences in fastq file
    with open(filepath, "r") as fastq_handle:
        for record in SeqIO.parse(fastq_handle, "fastq"):
            allseqs.append(str(record.seq))
    
#count instances of each barcode
    for entry in bc_seqs:
        rcomp = str(Seq(entry).reverse_complement())
        count = 0
        for line in allseqs:
            match = re.search(entry, line)
            match2 = re.search(rcomp, line)
            if match or match2:
                 count = count + 1
        bccounts.append((entry,count))
        
    
#dataframe of sequences and counts
    countsdf = pd.DataFrame(bccounts, columns=['seq', basename])
    
    readsum = countsdf[basename].sum()
    col1 = f"{basename}_share"
    countsdf[col1] = countsdf[basename]/readsum
    
    alldfs.append(countsdf)

#merge all timepoints
alldata = namesdf

for df in alldfs[0:]:
    alldata = alldata.merge(df, on="seq")

alldata.drop('seq', axis=1, inplace=True)


#Creates CSV of raw read counts
rawcounts = alldata.copy()

for col in rawcounts.columns:
    if col.endswith('_share'):
        rawcounts.drop(col, axis=1, inplace=True)

rawcounts.to_csv(f"bc_count_results/{rep}_rawcounts.csv", index=False)

#Creates CSV of read proportions (share of mapped reads)
props = pd.DataFrame()
for col in alldata.columns:
    if col == 'strain':
        props[col] = alldata[col]
    if col.endswith('_share'):
        colname = col.replace('_share', '')
        props[colname] = alldata[col]
        
props.to_csv(f'bc_count_results/{rep}_proportions.csv', index=False)

#Creates CSV of read proportions normalized to time0
normalize = props.copy()

time0 = normalize['0']

for col in normalize.columns:
    if col != 'strain':
        normalize[col] = normalize[col] / time0

normalize.to_csv(f'bc_count_results/{rep}_normalized.csv', index=False)

#Creates CSV of proportion of reads to sum of of S/D read counts for model input
toref = rawcounts.copy()

refdf = toref[toref['strain'].str.startswith(refstrain)].copy()
numcols = refdf.select_dtypes(include='number').columns
ref_sums = refdf[numcols].sum()
#drop ref rows
toref = toref[~toref['strain'].str.startswith('refstrain')]

#add ref barcode sums as a row
ref_row = {col: ref_sums[col] if col in numcols else ref if col == 'strain' else '' 
           for col in toref.columns}

toref = pd.concat([toref, pd.DataFrame([ref_row])], ignore_index=True)

#divide each timepoint by corresponding reference read sum (including reference sum itself)
toref[numcols] = toref[numcols].div(ref_sums)

toref.to_csv(f'bc_count_results/{rep}_prop_toref.csv', index=False)