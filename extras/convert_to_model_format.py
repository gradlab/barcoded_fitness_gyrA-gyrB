import os
import pandas as pd
import numpy as np

#Edit me
#Also edit hours/generations conversion (see line 49)
resultsdir = 'bc_count_results' #output of bc_count.py
genotypefile = 'genotypes_gyragyrb.txt' #barcodes and associated genotypes
filterend = 'prop_toref.csv' #replace with the type of results sheet you want (normalized, raw counts, etc.)
#filterstart = '2' #to get only ones from certain replicates - remove from master-sheet code if not necessary

#Make master-sheet with normalized values for each strain background in each competition replicate
masterlist = []

for file in os.listdir(resultsdir):
    if file.endswith(filterend): #and file.startswith(filterstart):
        basename, extension = os.path.splitext(file)
        rep = basename.split('_', 1)[0]
        file_path = os.path.join(resultsdir, file)
        df = pd.read_csv(file_path)
        df['replicate'] = rep 
        masterlist.append(df)

mastersheet = pd.concat(masterlist, ignore_index=True) 

#Add genotype
for index, row in mastersheet.iterrows():
    mastersheet['temp_parental'] = mastersheet['strain'].str.split('.').str[0]        
        
genotypedf = pd.read_csv(genotypefile, sep='\t')
genotypedict = genotypedf.set_index('strain')['genotype'].to_dict()
gyradict = genotypedf.set_index('strain')['gyra'].to_dict()
gyrbdict = genotypedf.set_index('strain')['gyrb'].to_dict()
mastersheet['genotype'] = mastersheet['temp_parental'].map(genotypedict)
mastersheet['gyra'] = mastersheet['temp_parental'].map(gyradict)
mastersheet['gyrb'] = mastersheet['temp_parental'].map(gyrbdict)
mastersheet.drop(['temp_parental'], axis=1, inplace=True)

#Create unique_id column
mastersheet['unique_id'] = mastersheet['strain'].astype(str) + '_' + mastersheet['replicate'].astype(str)

#Convert to long format
timepoints = [col for col in mastersheet.columns if col.isdigit()]

mastersheet_melt = pd.melt(mastersheet, id_vars=['strain','replicate','genotype','gyra','gyrb','unique_id'],
        value_vars=timepoints,
        var_name='time', value_name='prop_to_ref')

#Change from hours to generations - EDIT THIS!
    #(mastersheet_melt['time'] == TIMEPOINT)
    #choices = [CORRESPONDING # OF GENERATIONS IN REPLICATE ORDER]

mastersheet_melt['time'] = mastersheet_melt['time'].astype(float)

conditions = [
    (mastersheet_melt['replicate'] == 'A') & (mastersheet_melt['time'] == 2),
    (mastersheet_melt['replicate'] == 'B') & (mastersheet_melt['time'] == 2),
]
choices = [0.78, 1.42]

mastersheet_melt['time'] = np.select(conditions, choices, default=mastersheet_melt['time'])

conditions = [
    (mastersheet_melt['replicate'] == 'A') & (mastersheet_melt['time'] == 4),
    (mastersheet_melt['replicate'] == 'B') & (mastersheet_melt['time'] == 4),
]
choices = [3.77, 3.81]

mastersheet_melt['time'] = np.select(conditions, choices, default=mastersheet_melt['time'])

#Export
mastersheet_melt.rename(columns={'strain': 'barcode'}, inplace=True)
mastersheet_melt.to_csv('model_input.csv', index=False)