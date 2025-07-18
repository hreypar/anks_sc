#! /usr/bin/env python

# Definition of dataframes that will be used in other parts of the workflow
#
import pandas as pd

##### read in the SRA run table CSV file and store in dataframe
#
# use rename method on the dataframe to change the column names
sra_runtable = pd.read_csv(
    config['sra_runtable']['filename']
    ).rename(columns={
        config['sra_runtable']['sampid']: 'sampid',
	config['sra_runtable']['runid']: 'runid',
	config['sra_runtable']['layout']: 'layout'
    })

##### read in metadata table CSV file and store in dataframe
#
# maybe it wasnt necessary to change all these colnames
# but what if the names matched? we'd lose control of them
meta_table = pd.read_csv(
    config['samples_metadata']['filename']
    ).rename(columns={
        config['samples_metadata']['sampid_meta']: 'sampid_meta'
    })

##### Merge sra_runtable and meta_table based on sampid
#
# basically add metadata to sra_runtable
merged_sra_meta = sra_runtable.merge(meta_table, left_on='sampid', right_on='sampid_meta', how='left')
# I cant remember why lol

##### Create samples dataframe
#	
# group the rows of the sra_runtable dataframe by sampname_meta and apply the list function to runid column
# (because there could be more than one run by sample)
# the list function takes an iterable (here a column of the df) and returns a list containing all the elements of the iterable
# 
# this returns a new dataframe where for each sample, the runid column contains a list of all the SRA run IDs associated with that sample
samples = pd.DataFrame(merged_sra_meta.groupby('sample')['runid'].apply(list))


# set runid as the index of sra_runtable and store the resulting dataframe in runs 
# (runs is a new dataframe that has the same rows as sra_runtable but with runid column as the index)
# this means that we could use the .loc method to access rows/cols in the dataframe by index
runs = merged_sra_meta.set_index('runid')

# FINALLY make sample index for the meta table
meta_table.set_index('sample', inplace=True)

