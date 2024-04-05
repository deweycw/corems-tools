# Python script for aligning CoreMS outputs
# RMB Last updated  2/07/2024
# Contributors: Christian Dewey, Yuri Corilo, Will Kew,  Rene Boiteau

# Import the os module
import os
import pandas as pd
import time
import numpy as np
from coremstools import features as lcms





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
    data_dir='./testdata/'


    samplelist = [data_dir+f for f in os.listdir(data_dir) if '.csv' in f]

    features = lcms.Features()

    features.Align(samplelist)
    features.GapFill()
    featurelist.to_csv(data_dir + './feature_list.csv')
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
