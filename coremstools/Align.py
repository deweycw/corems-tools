from numpy import mean, std
from tqdm import tqdm
import numpy as np
import pandas as pd
import dask.dataframe as dd

from coremstools.Parameters import Settings

class Align:

    def Align_experimental(self, sample_list):
        """
        Method for assembling an aligned feature list. The aligned feature list is a dataframe containing a row for each [molecular formula]-[retention time] pair (what we call a feature) in the entire dataset. The dataframe contains the intensity of each feature in each sample in the data, as well as the average and stdev of each of the following parameters: measured m/z of the feature; calibrated m/z of the feature; resolving power of the instrument at the measured m/z; m/z error score; istopologue similarity score; confidence score; S/N; and dispersity. 

        Parameters 
        ----------
        sample_list : str
            Dataframe containing sample list. Must contain 'File' column with the name of each Thermo .raw file in the dataset. 
        """
        def ensure_same_columns(list_sample_csv):
            
            all_cols = []
            correct_order = []

            for sample_csv in list_sample_csv:
                
                print('\t' + sample_csv.split('/')[-1].split('.')[0])
                sample_cols = list(pd.read_csv(sample_csv).columns)
                if len(sample_cols) > len(correct_order):

                    correct_order = sample_cols

                temp = [col for col in sample_cols if col not in all_cols]
                
                all_cols = all_cols + temp
            
            temp_order = [c for c in correct_order if (c != 'Time') & (c != 'file') & ('Unnamed' not in c)]
            correct_order = ['file', 'Time'] + temp_order 

            for sample_csv in list_sample_csv:
                
                sample_temp = pd.read_csv(sample_csv)
                
                for col in all_cols:
                    
                    if col not in sample_temp.columns:
                        
                        sample_temp[col] = np.nan
                
                sample_temp = sample_temp[correct_order]
                sample_temp.to_csv(sample_csv, index=False)
            
        print('running alignment...')

        assignments_dir = Settings.assignments_directory
        
        list_sample_csv = [assignments_dir + f.replace('.raw', Settings.dispersity_addend + '.csv') for f in sample_list['File']]
        print('\tchecking columns...')
        ensure_same_columns(list_sample_csv)
        
        shared_columns = ['Time', 'Molecular Formula','Molecular Class', 'Ion Charge', 'Calculated m/z', 'Heteroatom Class',  'DBE']

        averaged_cols = ['m/z',
                    'm/z Error (ppm)',
                    'Calibrated m/z',
                    'Resolving Power',
                    'm/z Error Score',
                    'Isotopologue Similarity',
                    'Confidence Score',
                    'S/N',
                    'Dispersity']
        
        glob_str = assignments_dir + '*' + Settings.dispersity_addend + '.csv'

        all_results_read = dd.read_csv(list_sample_csv)

        all_results_shrink = all_results_read[all_results_read['Molecular Formula'].notnull()]
        
        def add_feature(row):

            z = row['Molecular Formula'] + str(row['Time'])
            return z
        
        all_results_shrink['feature'] = all_results_shrink.apply(add_feature, axis = 1) #['Molecular Formula'] + '--' + str(all_results_shrink['Time'])
        
        print('\tresetting index...')
        all_results = all_results_shrink.set_index('feature', sort = False)

        averaged_params = all_results[averaged_cols]

        averaged = averaged_params.groupby(by='feature').mean()

        stdev = averaged_params.groupby(by='feature').std()

        joined = averaged.join(stdev,lsuffix = '_mean', rsuffix = '_sd')

        shared = all_results[shared_columns]
        

        joined = joined.join(shared)

        for file in sample_list['File']:

            sub = all_results[all_results['file'] == file].copy()
            
            sub= sub.rename(columns = {'Peak Height':'Intensity in: ' + file})
             
            joined = joined.join(sub['Intensity in: ' + file])

        n_samples = all_results.groupby(by='feature').size()

        n_samples = n_samples.rename('N Samples')

        joined = joined.join(n_samples.to_frame(name='N Samples'))
        
        joined['feature'] = joined['Molecular Formula'] + str(joined['Time'])
        joined = joined.drop_duplicates(subset=['feature'])
        print('writing to .csv...')

        joined.to_csv(Settings.assignments_directory + 'feature_list.csv',index = False, single_file = True, header_first_partition_only = True)
        
        return joined
        

    def Align(self, sample_list):

        assignments_dir = Settings.assignments_directory
        disp_addend = Settings.dispersity_addend

        for file in sample_list['File']:

            print('  ' + file)

            file = assignments_dir + file.split('.')[0] + disp_addend + '.csv'

            results = dd.read_csv(file)
            
            results = results[results['Molecular Formula'].notnull()]
            
            #results['feature'] = list(zip(results['Time'],results['Molecular Formula']))
            results['feature'] = str(results['Time']) + results['Molecular Formula']
            
            file_name = file.replace('.csv','').split('/')[-1]

            masterresults['Intensity: '+file_name]={}
            
            pbar = tqdm(range(len(results)), ncols=100)

            for ix in pbar:

                row = results.iloc[ix,:]
                
                if row['feature'] not in masterresults['Time'].keys():

                    for col in shared_columns:
                        
                        masterresults[col][row['feature']] = row[col]

                    current_elements = [x.rstrip('0123456789') for x in row['Molecular Formula'].split()]
                    
                    for element in current_elements:

                        if element not in elements:

                            elements.append(element)

                            masterresults[element]={}

                        masterresults[element][row['feature']]=row[element]

                    masterresults['Intensity: ' + file_name][row['feature']] = int(row['Peak Height'])

                    for c in averaged_cols:

                        masterresults[c][row['feature']] = [row[c]]

                else:
                    masterresults['Intensity: ' + file_name][row['feature']] = int(row['Peak Height'])
                    
                    for c in averaged_cols:
                        
                        masterresults[c][row.feature].append(row[c])
                
        print('  writing N Samples column')

        pbar = tqdm(masterresults['m/z'].keys(), ncols=100)

        for key in pbar:

            masterresults['N Samples'][key] = len(masterresults['m/z'][key])

            for c in averaged_cols:

                masterresults[c+' stdev'][key] = std(masterresults[c][key])
                masterresults[c][key] = mean(masterresults[c][key])
        results_df = DataFrame(masterresults).fillna(0)
        cols_at_end = [c for c in results_df.columns if 'Intensity' in c ]
        results_df = results_df[[c for c in results_df if c not in cols_at_end] + [c for c in cols_at_end if c in results_df]]
        
        return(results_df)
        

