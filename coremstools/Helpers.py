import pandas as pd
import numpy as np
import tqdm
import warnings

warnings.filterwarnings("ignore")
import sys
sys.path.append("./")

"""
Helper functions for CoreMS assignments.

Parameters 
----------
sample_df : DataFrame 
    Pandas DataFrame containing a 'File' column with the name of each .raw file in the dataset (not full path). Defaults to None.
t_interval : float
    Interval (min) over which scans are averaged in CoreMS LC assignments. Defaults to 2 (min).   

Methods
-------
add_mol_class()
    Adds molecular class to CoreMS assignments & creates a 'Molecular Class' column. Altered assignment .csvs are rewritten with new columns. 
run_internal_std_qc(timerange=[10,12]) -> DataFrame
    Runs the quality control checks with internal standard m/z. Returns copy of sample list DataFrame with additional columns for internal standard area, retention time, and QC pass/fail flag.
run_assignment_error_plot()
    For each sample in the sample list, this method creates .jpg plots of (i) m/z Error (ppm) v. m/z and (ii) Molecular Classes of assignments over separation. The .jpgs are saved in the directory defined by Settings.assignments_directory.
run_dispersity_calculation()
    Runs dispersity calculation on each m/z in the CoreMS assignment file corresponding to each sample. The CoreMS assignment files are copied and saved as [SAMPLE_NAME] + dispersity_addend +'.csv' in the directory defined by Settings.assignments_directory. ***** Currently quite slow. Would be good to do this calculation after the feature list is assembled. ******
"""
def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def get_heteroatoms(df):
    cols = df.columns
    add = False
    heteros = []
    for c in cols:
        if add:
            if has_numbers(c):
                continue
            elif c == 'Time':
                break
            else:
                heteros.append(c)
        if c == 'O':
            add = True
        elif c == 'Time':
            add = False
    return heteros

def assign_mol_class(complete_results_df,molclasses):
                # creates list of mol_classes based on heteroatom list
    all_results = complete_results_df    # adds 'mol_class' column

    i = 0 # iterable 
    p = 0 # max len holder
    j = 0 # index of max len mol class
    
    for i in range(0,len(molclasses)):

        mc = molclasses[i]
        
        if (len(mc) > p) and (mc != 'Unassigned'):
            
            p = len(mc)
            
            j = i

    all_elements = get_elements(molclasses[j])

    all_results['ID'] = range(0,len(all_results))

    times = all_results['Time'].unique()

    holder = []
    sizenp = 0
    pbar = tqdm.tqdm(times)

    for t in pbar:

        time_average = all_results[all_results['Time'] == t]

        sizenp = sizenp + len(time_average)

        for m in molclasses:

            if m != 'Unassigned':
                elements = get_elements(m)

                sub = get_molclass_subset(elements, all_elements,time_average[~time_average['Molecular Formula'].isna()]) 
                sub['Molecular Class'] = m    


            elif m == 'Unassigned':
                sub = time_average[time_average['Molecular Formula'].isna()] 
                sub['Molecular Class'] = m

            #print('\t%s: %s' %(m,len(sub)))

            pbar.set_description_str(desc="Assigning molecular class %s at time %s" % (m,t) , refresh=True)

            holder.append(sub)


    results = pd.concat(holder)
    return results 


def _get_mol_classes(add, base):
    
    new = []
    remain = []
    
    new.append(base)

    for i in range(len(add)):
        
        new.append(base + add[i])

        new2 = []

        remain = add[i+1:]

        for j in remain:

            new2 = add[i] + j

            new.append(base + new2)

    return(new)


def get_mol_class(add):
    base = 'CHO'
    molclasses = []
    for i in range(len(add)):
        e = add[i]
        temp = _get_mol_classes(add[i:], base = base)
        base = base+e
        molclasses = molclasses + temp


    output = []
    for x in molclasses:
        if x not in output:
            output.append(x)

    output.append('Unassigned')
    return output


def get_elements(molclass):

    import re
    
    elements = [] 

    elements = re.findall('[A-Z][^A-Z]*', molclass)
    
    return elements


def get_molclass_subset(included_elements, all_elements, all_results):
    
    tdf = all_results

    for e in all_elements:

        try:
        
            tdf[e].fillna(0, inplace = True)

        except:

            pass

    excluded_elements = [e for e in all_elements if e not in included_elements]
    
    for e in included_elements:
        
        try:

            tdf = tdf[tdf[e]>0]

            for j in excluded_elements:

                try:
                    
                    tdf = tdf[tdf[j]==0]

                except:

                    pass

        except:
            
            pass

    return tdf


def get_ratios(results):

    results[results['Is Isotopologue']==0]['O/C'] = results['O'] / results['C']
    results[results['Is Isotopologue']==0]['H/C'] = results['H'] / results['C']
    results[results['Is Isotopologue']==0]['N/C'] = results['N'] / results['C']

    return results 


def add_mzwindow_col(df):    

    df['m/z window'] = df.index
    df['Window Size (m/z)'] = df.index
    pbar = tqdm.tqdm(zip(df['file'], range(len(df['file']))))

    for file, r in pbar:

        if 'StdMix' not in file:

            if ('400_500' in file) or ('400-500' in file):

                df['m/z window'].iloc[r] = '400-500 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"

            elif ('500_600' in file) or ('500-600' in file):

                df['m/z window'].iloc[r] = '500-600 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"
        
            elif ('600_700' in file) or ('600-700' in file):

                df['m/z window'].iloc[r] = '600-700 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"

            elif ('700_800' in file) or ('700-800' in file):

                df['m/z window'].iloc[r] = '700-800 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"
            
            elif ('300_500' in file) or ('300-500' in file):

                df['m/z window'].iloc[r] = '300-500 m/z'
                df['Window Size (m/z)'].iloc[r] = "200"

            
            elif ('400_600' in file) or ('400-600' in file):

                df['m/z window'].iloc[r] = '400-600 m/z'
                df['Window Size (m/z)'].iloc[r] = "200"
            
            elif ('600_800' in file) or ('600-800' in file):

                df['m/z window'].iloc[r] = '600-800 m/z'
                df['Window Size (m/z)'].iloc[r] = "200"
                
            elif 'full' in file:

                df['m/z window'].iloc[r] = '200-1200 m/z'
                df['Window Size (m/z)'].iloc[r] = "1000"
        
        else:
            df['m/z window'].iloc[r] = '200-1200 m/z'
            df['Window Size (m/z)'].iloc[r] = "1000"
        
        pbar.set_description_str(desc="Assigning window size columns for %s" %(file.split('/')[-1]) , refresh=True)

    return df 


def getUniqueFeatures(df):    
    '''
    Notes by Christian Dewey, 23-Feb-27

    OVERVIEW:
    
        - subset m/z assignments by time bin
        - sort subset by m/z error (ascending)
        - remove duplicate m.f. assignments, preserving the assignment with lowest m/z error (first hit on sorted subset), save as 'currentunique'
        - set the index of 'currentunique' to the 'Molecular Formula' column
        - for each file in the raw file list:
            - save subset the sorted m/z assignment list (with all files and duplicates) by file name as 'current_file'
            - remove duplicate formulae from file subset copies
            - rename 'Peak Height' to filename for each file
            - set 'Molecular Formula' to index for copy with renamed 'Peak Height' col
            - join the renamed 'Peak Height' col to 'currentunique'

    RETURNS:    Pandas dataframe 
                Contains unique m.f. assignments in each time bin, and intenisty of unique m.f. feature in each file within each time bin

    
    '''
    uniquelist=[]
    for time in df.Time.unique():
        current=df[df.Time==time]
        current=current.sort_values(by=['m/z Error (ppm)'],ascending=True)
        currentunique=current.drop_duplicates(subset=['Molecular Formula'])
        currentunique=currentunique[currentunique['C']>1]
        currentunique=currentunique.set_index(['Molecular Formula'],drop=False)
        for file in df['file'].unique():
            current_file=current[current['file']==file].drop_duplicates(subset=['Molecular Formula'])
            current_file=current_file.rename(columns={'Peak Height':file})
            current_file=current_file.set_index(['Molecular Formula'],drop=False)
            currentunique=currentunique.join(current_file[file])
        '''for mzw in df['Window Size (m/z)'].unique():
            current_file=current[current['Window Size (m/z)']==mzw].drop_duplicates(subset=['Molecular Formula'])
            wlbl = mzw + ' m/z window'
            current_file=current_file.rename(columns={'Peak Height':wlbl})
            current_file=current_file.set_index(['Molecular Formula'],drop=False)
            currentunique=currentunique.join(current_file[wlbl])'''
        uniquelist.append(currentunique)

    unique_results=pd.concat(uniquelist,ignore_index=True)
    unique_results['N/C']=unique_results['N']/unique_results['C']

    return unique_results




def addRepCol(data_df):

    data_df['Rep'] = data_df.index


    for file in data_df['file'].unique():

        #print(file)

        if ('rep2' in file) or ('_02.' in file):

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 2
            data_df[data_df['file'] == file] = temp


        else:

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 1
            data_df[data_df['file'] == file] = temp

    #print(data_df['Rep'].unique())
    return data_df 


def add_mz_window_colsl(data_df):

    data_df['Rep'] = data_df.index

    pbar = tqdm.tqdm(data_df['file'].unique())
           
    for file in pbar:

        if ('rep2' in file) or ('_02.' in file):

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 2
            data_df[data_df['file'] == file] = temp


        else:

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 1
            data_df[data_df['file'] == file] = temp
        pbar.set_description_str(desc="Adding rep column to %s" %(file) , refresh=True)
    return data_df 


def blankSubtract(df, blnkthresh = 0.8):
    # must be performed on df with unique assignments 
    holder = []
    for file in df['file'].unique():
        
        sub = df[df['file'] == file]

        blkf = sub['blank file'].iloc[0]

        sub[sub[file]== np.nan] = 0  # each file column contains intensities of feature; if feature is not present, nan assigned; need to convert to zero for blank subtract

        nom = sub[file]
        den = sub[blkf]

        nom = nom.replace(np.nan,0)  # features not present in sample
        den = den.replace(np.nan,1)  # features not present in blank

        if file != blkf:
            nom = nom
        elif file == blkf:
            nom = nom * (blnkthresh * 0.8)  # multiplication enables removal of these features from blank file 

        sub['blank subtract'] = nom/den  # ratio of intensities of features in sample and blank

        holder.append(sub)

    df_end = pd.concat(holder)

    df_end = df_end[df_end['blank subtract'] > blnkthresh]  # only select features that do not appear in blanks

    return df_end


def repCombine(df):

    for file in df['file'].unique():

        df[df[file] == np.nan] = 0

        if 'rep2' not in file:

            if '.raw' in file:

                rep2file = file.split('.')[0]+'_rep2.raw'

            else:

                rep2file = file + '_rep2'
            
            avfile = file + '_av'
            
            df[avfile] = (df[file] + df[rep2file]) / 2

    return df


def normMS(df,fulldf):

    max_i = max(fulldf['Peak Height'].values)

    df['Normalized Peak Height'] = df['Peak Height'] / max_i

    return df

