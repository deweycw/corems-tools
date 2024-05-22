from numpy import mean, std
from tqdm import tqdm
import numpy as np
import dask.dataframe as dd
from coremstools.Parameters import Settings

class Align:

    def Align(self, sample_list):
        """
        Method for assembling an aligned feature list. The aligned feature list is a dataframe containing a row for each [molecular formula]-[retention time] pair (what we call a feature) in the entire dataset. The dataframe contains the intensity of each feature in each sample in the data, as well as the average and stdev of each of the following parameters: measured m/z of the feature; calibrated m/z of the feature; resolving power of the instrument at the measured m/z; m/z error score; istopologue similarity score; confidence score; S/N; and dispersity. 

        Parameters 
        ----------
        sample_list : str
            Dataframe containing sample list. Must contain 'File' column with the name of each Thermo .raw file in the dataset. 
        """

        assignments_dir = Settings.assignments_directory

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
        
        print('running alignment ...')
                
        glob_str = assignments_dir + '*' + Settings.dispersity_addend + '.csv'

        all_results = dd.read_csv(glob_str)

           
        all_results = all_results[all_results['Molecular Formula'].notnull()]

        all_results['feature'] = all_results['Molecular Formula'] + '--' + str(all_results['Time'])

        all_results = all_results.set_index('feature')

        #files = list(sample_list['File'])

        #features = all_results.drop_duplicates(subset = ['feature'])['feature']

        #features.set_index('feature', inplace = True)

        averaged_params = all_results[averaged_cols].copy()

        averaged = averaged_params.groupby(by='feature').mean()

        stdev = averaged_params.groupby(by='feature').std()

        joined = averaged.join(stdev,lsuffix = '_mean', rsuffix = '_sd')

        shared = all_results[shared_columns].copy()

        joined = joined.join(shared)

        for file in sample_list['File']:

            print(file)
            
            sub = all_results[all_results['file'] == file].copy()
            
            sub = sub.rename(columns = {'Peak Height':'Intensity in: ' + file})
            
            joined = joined.join(sub['Intensity in: ' + file])

        n_samples = all_results.groupby(by='feature').size()

        n_samples = n_samples.rename('N Samples')

        joined = joined.join(n_samples.to_frame(name='N Samples'))

        print(joined.columns)

        print('writing to .csv...')

        joined.to_csv(Settings.assignments_directory + 'feature_list.csv', single_file = True)

        return joined
        
    def old():            
        
        pass 
        '''print('  writing N Samples column')

        pbar = tqdm(masterresults['m/z'].keys(), ncols=100)

        for key in pbar:

            masterresults['N Samples'][key] = len(masterresults['m/z'][key])

            for c in averaged_cols:

                masterresults[c+' stdev'][key] = std(masterresults[c][key])
                masterresults[c][key] = mean(masterresults[c][key])
        results_df = DataFrame(masterresults).fillna(0)
        cols_at_end = [c for c in results_df.columns if 'Intensity' in c ]
        results_df = results_df[[c for c in results_df if c not in cols_at_end] + [c for c in cols_at_end if c in results_df]]
        
        return(results_df)'''
    


def old(self, sample_list):

        assignments_dir = Settings.assignments_directory
        disp_addend = Settings.dispersity_addend

        for file in sample_list['File']:

            print('  ' + file)

            file = assignments_dir + file.split('.')[0] + disp_addend + '.csv'

            results = dd.read_csv(file)
            
            results = results[results['Molecular Formula'].notnull()]
            
            results['feature'] = list(zip(results['Time'],results['Molecular Formula']))
            
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
        

