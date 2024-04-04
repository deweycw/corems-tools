# Python script for aligning CoreMS outputs
# RMB Last updated  2/07/2024
# Contributors: Christian Dewey, Yuri Corilo, Will Kew,  Rene Boiteau

# Import the os module
import os
import pandas as pd
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from coremstools import features as lcms




def gapfill(featurelist):
    featurelist['gapfill']=False
    featurelist['gapfill flag']=False
    for i, row in featurelist.iterrows():
        resolution=row['Resolving Power']
        mass=row['Calibrated m/z']
        time=row['Time']
        mrange=[mass*(1-2/resolution),mass*(1+2/resolution)]
        matches=featurelist[(featurelist['Calibrated m/z']>mrange[0])&(featurelist['Calibrated m/z']<mrange[1])&(featurelist['Time']==time)]
        if(len(matches)>1):
            featurelist.loc[i,'gapfill']=True
            featurelist.loc[i,featurelist.filter(regex='Intensity').columns]=matches.filter(regex='Intensity').sum(axis=0)
            if featurelist.loc[i,'Confidence Score']<max(matches['Confidence Score']):
                featurelist.loc[i,'gapfill flag']=True
    return(featurelist)

def gapfill2(featurelist):
    """
    Fills gaps in 'featurelist' pandas DataFrame intensities based on calibrated m/z and resolving power.

    Args:
        featurelist: A pandas DataFrame with columns 'Resolving Power', 'Calibrated m/z',
                    'Time', and multiple 'Intensity' columns.

    Returns:
        A pandas DataFrame with a new 'gapfill' column indicating if a gap was filled
        and the filled intensity values (sum of matching rows) in the 'Intensity' columns.
    """

    resolution = featurelist['Resolving Power']
    mass = featurelist['Calibrated m/z']
    time = featurelist['Time']
    mrange = [(1 - 2/resolution) * mass, (1 + 2/resolution) * mass]

    # Vectorized filtering: find rows within mass range and matching time
    matches = featurelist[(featurelist['Calibrated m/z'] >= mrange[0]) &
                        (featurelist['Calibrated m/z'] <= mrange[1]) &
                        (featurelist['Time'] == time)]

    print(matches.columns)

    # Vectorized gap filling: create new DataFrame with gapfilled intensities
    intensity_cols = featurelist.filter(like='Intensity').columns
    gapfilled_df = (matches.groupby(['Calibrated m/z', 'Time'])
                    .agg(gapfill=(matches[intensity_cols] == 0).any(), *[(col, 'sum') for col in intensity_cols])
                    .reset_index()
    )

   # gapfilled_df = (

   #     matches.groupby(['Calibrated m/z', 'Time']).agg(gapfill=(matches[intensity_cols] == 0).any(), *[(col, 'sum') for col in intensity_cols])
   #     .reset_index()
   # )


    # Update original DataFrame with gapfill flag and summed intensities
    #featurelist.update(gapfilled_df[['gapfill']])
    featurelist.update(gapfilled_df[['gapfill']]
                      .merge(gapfilled_df[intensity_cols].groupby(['Resolving Power', 'Calibrated m/z', 'Time']).sum(),
                            on=['Resolving Power', 'Calibrated m/z', 'Time'],
                            how='left'))

    return featurelist

def OHNratios(featurelist):
    # Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
    featurelist['O/C']=featurelist['O']/featurelist['C']
    featurelist['H/C']=featurelist['H']/featurelist['C']
    featurelist['N/C']=featurelist['N']/featurelist['C']
    featurelist['NOSC'] =  4 -(4*featurelist['C']
                            + featurelist['H']
                            - 3*featurelist['N']
                            - 2*featurelist['O'])/featurelist['C']

def mz_error_flag(featurelist):
    featurelist=featurelist.sort_values(by=['Calculated m/z'])
    featurelist['rolling error']=featurelist['m/z Error (ppm)'].rolling(int(len(featurelist)/50),center=True,min_periods=0).mean()
    featurelist['mz error flag']=abs(featurelist['rolling error']-featurelist['m/z Error (ppm)'])/(4*featurelist['m/z Error (ppm) stdev'])
    return(featurelist)

def blank_flag(featurelist):
    featurelist['Max Intense']=featurelist.filter(regex='Intensity').max(axis=1)
    featurelist['blank']=featurelist['Intensity:'+blankfile.replace(dfiletype,'')].fillna(0)/featurelist['Max Intense']
    return(featurelist)

if __name__ == '__main__':

    #### Change file settings here
    global data_dir
    data_dir='/Users/christiandewey/Code/coremstools/test/testdata/'

    global dfiletype
    dfiletype='.raw'

    ##### End user input

    starttime = time.time()

    #featurelist=pd.read_csv(data_dir+featurelist_file)

    #tracemalloc.start()

    samplelist = [data_dir+f for f in os.listdir(data_dir) if '.csv' in f]

    features = lcms.Features()

    features.Align(samplelist)

    print(features.results)
    '''featurelist=featurelist_aligner(samplelist)
    #featurelist.to_csv(data_dir + 'feature_list.csv')
    featurelist=gapfill(featurelist)
    featurelist.to_csv(data_dir + 'feature_list.csv')'''
    #OHNratios(featurelist)
    '''

    print("# Features in list: " + str(len(featurelist)))

    ### Remove features detected in the blank within 50% of the max intensity.

    featurelist=blank_flag(featurelist)
    print("# Features, blank corrected: " + str(len(featurelist[featurelist['blank']<0.5])))
    #featurelist=featurelist[featurelist['blank']<0.5]

    ### Remove features with average m/z more than 4x standard deviation of mean error.
    featurelist=mz_error_flag(featurelist)
    print("Unique results, error corrected: " + str(len(featurelist[featurelist['mz error flag']<1])))
    #featurelist=featurelist[featurelist['mz error flag']==True]

    print("Unique molecular formula: " + str(len(featurelist['Molecular Formula'].unique())))

    featurelist.to_csv(data_dir+featurelist_file)

    #print('\n\nCurrent memory in use: %.2f MB\nMaximum memory used: %.2f MB' %(tracemalloc.get_traced_memory()[0]/1000/1000,tracemalloc.get_traced_memory()[1]/1000/1000))
    #tracemalloc.stop()

    print('Total execution time: %.2f min' %((time.time()-starttime)/60))
    '''
