from numpy import mean, std
from tqdm import tqdm
import numpy as np
import pandas as pd
#import dask.dataframe as dd
#from dask.diagnostics import ProgressBar
#ProgressBar().register()

from coremstools.Parameters import Settings


class Align:

    def run(self, sample_list, include_dispersity = True):
        """
        Method for assembling an aligned feature list. The aligned feature list is a dataframe containing a row for each [molecular formula]-[retention time] pair (what we call a feature) in the entire dataset. The dataframe contains the intensity of each feature in each sample in the data, as well as the average and stdev of each of the following parameters: measured m/z of the feature; calibrated m/z of the feature; resolving power of the instrument at the measured m/z; m/z error score; istopologue similarity score; confidence score; S/N; and dispersity. 

        Parameters 
        ----------
        sample_list : str
            Dataframe containing sample list. Must contain 'File' column with the name of each Thermo .raw file in the dataset. 
        """
        def build_masterresults_dict(shared_columns, averaged_cols):
            
            masterresults={}

            for col in shared_columns:
                
                masterresults[col] = {}
            
            masterresults['N Samples'] = {}
            
            for col in averaged_cols:
                
                masterresults[col] = {}
                
                masterresults[col + '_se'] = {}

            return masterresults

        

        assignments_dir = Settings.assignments_directory

        shared_columns = ['Time','Molecular Formula',  'Calculated m/z', 'DBE', 'Is Isotopologue', 'Molecular Class' ,'Heteroatom Class']

        averaged_cols = ['m/z',
                    'm/z Error (ppm)',
                    'Calibrated m/z',
                    'Resolving Power',
                    'Confidence Score',
                    'm/z Error Score',
                    'Isotopologue Similarity',
                    'S/N',
                    'Dispersity',
                    'Retention Time']
        
        if include_dispersity is False:

            averaged_cols = ['m/z',
                    'm/z Error (ppm)',
                    'Calibrated m/z',
                    'Resolving Power',
                    'Confidence Score',
                    'm/z Error Score',
                    'Isotopologue Similarity',
                    'S/N']

        print('running alignment on ...')
                
        elements=[]

        masterresults = build_masterresults_dict(shared_columns, averaged_cols)
        used_elements = []

        for file in sample_list['File']:

            print('  ' + file)

            file = assignments_dir + file.split('.')[0] + '.csv'

            results = pd.read_csv(file)
            
            results = results[results['Molecular Formula'].notnull()]
            
            results['feature'] = list(zip(results['Time'],results['Molecular Formula']))
            
            file_name = file.replace('.csv','').split('/')[-1]

            masterresults['Intensity: '+file_name]={}
            
            pbar = tqdm(range(len(results)))

            for ix in pbar:

                row = results.iloc[ix,:]
                
                if row['feature'] not in masterresults['Time'].keys():

                    for col in shared_columns:
                        
                        masterresults[col][row['feature']] = row[col]

                    current_elements = [x.rstrip('0123456789') for x in row['Molecular Formula'].split()]
                    
                    for element in current_elements:

                        if element not in elements:

                            elements.append(element)
                            used_elements.append(element)

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

        pbar = tqdm(masterresults['m/z'].keys())

        for key in pbar:

            masterresults['N Samples'][key] = len(masterresults['m/z'][key])

            for c in averaged_cols:

                masterresults[c+'_se'][key] = std(masterresults[c][key]) / masterresults['N Samples'][key]
                masterresults[c][key] = mean(masterresults[c][key])

        results_df = pd.DataFrame(masterresults).fillna(0)
        cols_at_end = [c for c in results_df.columns if 'Intensity' in c ]
        
        final_col_list = shared_columns + [ f for f in averaged_cols] + [ f + '_se' for f in averaged_cols] 

        final_col_list = [f for f in final_col_list if (f != 'file') & (f != 'Peak Height')] + used_elements + ['N Samples'] + cols_at_end
        
        results_df = results_df[final_col_list]
        
        return(results_df)


    def Align_exp(self, sample_list, include_dispersity = True):
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
            rewrite = False
            first = True 
            for sample_csv in list_sample_csv:
                sample_cols = list(pd.read_csv(sample_csv).columns)
                if len(sample_cols) > len(correct_order):

                    correct_order = sample_cols

                temp = [col for col in sample_cols if col not in all_cols]

                if (len(temp) > 0) & (not first):
                    rewrite = True      
                
                first = False
                all_cols = all_cols + temp
            
            temp_order = [c for c in correct_order if (c != 'Time') & (c != 'file') & ('Unnamed' not in c)]
            correct_order = ['file', 'Time'] + temp_order 

            if rewrite:
                #print('\trewriting with updated order')
                for sample_csv in list_sample_csv:
                    
                    sample_temp = pd.read_csv(sample_csv)
                    
                    for col in all_cols:
                        
                        if col not in sample_temp.columns:
                            
                            sample_temp[col] = np.nan
                    
                    sample_temp = sample_temp[correct_order]
                    sample_temp.to_csv(sample_csv, index=False)
            
        #print('running alignment...')

        assignments_dir = Settings.assignments_directory
        
        list_sample_csv = [assignments_dir + f.replace('.raw', '.csv') for f in sample_list['File']]
        #print('checking columns...')
        ensure_same_columns(list_sample_csv)
        
        #shared_columns = ['Time', 'Molecular Formula','Molecular Class', 'Ion Charge', 'Calculated m/z', 'Heteroatom Class',  'DBE']
        
        shared_columns = ['file','Peak Height','Time', 'Molecular Formula', 'Ion Charge', 'Calculated m/z', 'DBE']

        averaged_cols = ['m/z',
                         'm/z Error (ppm)',
                         'Calibrated m/z',
                         'Resolving Power',
                         'Confidence Score',
                         'S/N',
                         'Dispersity']
        
        glob_str = assignments_dir + '*' +  '.csv'

        all_results_read = dd.read_csv(list_sample_csv)

        all_results_shrink = all_results_read[all_results_read['Molecular Formula'].notnull()]
        
        def add_feature(row):

            z = str(row['Time']) + '--' + row['Molecular Formula'] 
            return z
        
        all_results_shrink['feature'] = all_results_shrink.apply(add_feature, axis = 1) #['Molecular Formula'] + '--' + str(all_results_shrink['Time'])
        
        #print('resetting index...')
        all_results = all_results_shrink.set_index('feature', sort = False)

        averaged_params = all_results[averaged_cols]

        averaged = averaged_params.groupby(by='feature').mean()

        stdev = averaged_params.groupby(by='feature').std()

        joined = averaged.join(stdev,lsuffix = '_mean', rsuffix = '_se')

        shared = all_results[shared_columns]
        
        joined = joined.join(shared)
              
        #print('assembling intensities...')
        flist = list(sample_list['File'])
        def assemble_intensities(group):
            file_keys = [f.split('/')[-1] for f in list(group['file'])]
            peak_heights = list(group['Peak Height'])
            missing = [f for f in flist if f not in file_keys]
            add_dict = {m:0 for m in missing}
            int_dict = {k : int(i) for k, i in zip(file_keys, peak_heights)}
            full_dict = {**int_dict, **add_dict}
            for f in flist:
                group[f] = full_dict[f]
            return group
        
        feature_groups = joined.groupby('feature').apply(assemble_intensities)
        
        n_samples = feature_groups.groupby(by='feature').size()
        n_samples = n_samples.rename('N Samples')

        joined2 = feature_groups.join(n_samples.to_frame(name='N Samples'))
        joined3 = joined2.groupby(joined2.index).first() #  [last(joined2.index.drop_duplicates())]

        final_col_list = [ f + '_mean' for f in averaged_cols] + [ f + '_se' for f in averaged_cols] + ['N Samples']
        final_col_list = shared_columns + final_col_list
        final_col_list = [f for f in final_col_list if (f != 'file') & (f != 'Peak Height')] + flist
        return joined3[final_col_list]
        