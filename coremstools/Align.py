from numpy import mean, std, sqrt, nan
from tqdm import tqdm
import numpy as np
import pandas as pd

from coremstools.Parameters import Settings

import re 

class Align:

    def run(self, sample_list):
        """
        Method for assembling an aligned feature list. The aligned feature list is a dataframe containing a row for each [molecular formula]-[retention time] pair (what we call a feature) in the entire dataset. The dataframe contains the intensity of each feature in each sample in the data, as well as the average and stdev of each of the following parameters: measured m/z of the feature; calibrated m/z of the feature; resolving power of the instrument at the measured m/z; m/z error score; istopologue similarity score; confidence score; S/N; and dispersity. 

        Parameters 
        ----------
        sample_list : str
            Dataframe containing sample list. Must contain 'File' column with the name of each Thermo .raw file in the dataset. 
        """
        assignments_dir = Settings.assignments_directory

        shared_columns_def = ['Time', 'Molecular Formula', 'Calculated m/z', 'DBE', 'Is Isotopologue', 'Molecular Class', 'Heteroatom Class']
        
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
        
        all_dfs = []
        
        print('Reading and preprocessing files...')
        for file_path_suffix in tqdm(sample_list['File'], desc="Processing files"):
            file_name_no_ext = file_path_suffix.split('.')[0]
            full_file_path = f"{assignments_dir}{file_name_no_ext}.csv"

            try:
                df = pd.read_csv(full_file_path)
            except FileNotFoundError:
                print(f"Warning: File not found {full_file_path}, skipping.")
                continue

            df = df[df['Molecular Formula'].notnull()].copy()
            if df.empty:
                continue

            df['feature'] = df['Time'].astype(str) + '--' + df['Molecular Formula']
            df['file_name_for_pivot'] = 'Intensity: ' + file_name_no_ext
            all_dfs.append(df)

        if not all_dfs:
            print("No data to process after reading files.")
            return pd.DataFrame()

        print('Concatenating data...')
        concatenated_df = pd.concat(all_dfs, ignore_index=True)
        
        # Determine actual element columns based on formulas and presence in DataFrame
        unique_mf_series = concatenated_df['Molecular Formula'].drop_duplicates()
        parsed_elements_from_formulas = set()
        element_parser_re = re.compile(r'(\d*[A-Z][a-z]?)') # Parses elements including isotopes (e.g., C, 13C, Cl, 35Cl)
        for mf_str in unique_mf_series:
            if pd.notna(mf_str):
                parsed_elements_from_formulas.update(element_parser_re.findall(mf_str))
        
        actual_element_cols = sorted([el for el in parsed_elements_from_formulas if el in concatenated_df.columns])
        
        print('Grouping and aggregating data...')
        grouped = concatenated_df.groupby('feature')
        agg_funcs = {}

        #First, create columns for shared columns and actual element columns
        for col in shared_columns_def:
            if col in concatenated_df.columns: agg_funcs[col] = 'first'
        for col in actual_element_cols:
            if col in concatenated_df.columns: agg_funcs[col] = 'first'

        # Now, handle averaged columns
        averaged_cols_def = [col for col in averaged_cols if col in concatenated_df.columns]

        for col in averaged_cols_def:
            if col in concatenated_df.columns:
                # Only add 'mean' and 'std' aggregations here
                agg_funcs[col] = ['mean', 'std']

        if not agg_funcs and 'feature' in concatenated_df:
            aggregated_df = pd.DataFrame(index=concatenated_df['feature'].unique())
        elif not agg_funcs:
            print("No data to aggregate.")
            return pd.DataFrame()
        else:
            aggregated_df = grouped.agg(agg_funcs)
            aggregated_df.columns = ['_'.join(col).strip('_') for col in aggregated_df.columns.values]

        # Create 'N Samples' column using grouped.size()
        aggregated_df['N Samples'] = grouped.size()

        # Create Standard error columns

        for col in averaged_cols_def:
            std_col = col + '_std'
            se_col = col + '_se'
            if std_col in aggregated_df.columns and 'N Samples' in aggregated_df.columns:
                n_samples = aggregated_df['N Samples']
                aggregated_df[se_col] = aggregated_df[std_col] / sqrt(n_samples)

        rename_map = {}

        for col in shared_columns_def:
            processed_name = col + '_first'
            if processed_name in aggregated_df.columns:
                rename_map[processed_name] = col
        for col in actual_element_cols:
            processed_name = col + '_first'
            if processed_name in aggregated_df.columns:
                rename_map[processed_name] = col
        for col in averaged_cols_def:
            processed_name = col + '_mean'
            if processed_name in aggregated_df.columns:
                rename_map[processed_name] = col

        aggregated_df.rename(columns=rename_map, inplace=True)

        print('Getting intensities...')
        intensity_df = concatenated_df.pivot_table(index='feature', 
                                                   columns='file_name_for_pivot', 
                                                   values='Peak Height', 
                                                   fill_value=0)

        print('Joining aggregated data with intensities...')

        final_df = aggregated_df.join(intensity_df.astype(int), how='left')

        # Reorder columns to have shared columns, averaged columns, element columns, and then intensities
        ordered_cols_list = []
        for col in shared_columns_def:
            if col in aggregated_df.columns: # Check if the renamed column exists in aggregated_df
                ordered_cols_list.append(col)
        
        
        # 2. Averaged columns (now renamed mean, std, se)
        for col in averaged_cols_def:
            mean_col = col
            std_col = col + '_std'
            se_col = col + '_se'
            if mean_col in aggregated_df.columns: ordered_cols_list.append(mean_col)
            #if std_col in aggregated_df.columns: ordered_cols_list.append(std_col)
            if se_col in aggregated_df.columns: ordered_cols_list.append(se_col)

        # 3. Element columns (now renamed)
        for col in actual_element_cols:
            if col in aggregated_df.columns: # Check if the renamed column exists in aggregated_df
                ordered_cols_list.append(col)

        # 4. N Samples and Intensity columns
        if 'N Samples' in aggregated_df.columns: ordered_cols_list.append('N Samples')

        ordered_cols_list.extend(sorted(intensity_df.columns.tolist()))

        # Ensure the final DataFrame only contains the columns in the ordered list
        ordered_cols_list_present = [col for col in ordered_cols_list if col in final_df.columns]
        # Reindex final_df to the desired order
        final_df = final_df[ordered_cols_list_present]

        # Fill NaN values with 0
        final_df = final_df.fillna(0)

        print('Alignment finished.')
        return final_df
    
    def run_legacy(self, sample_list, include_dispersity = True):
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
                    'S/N',
                    'Retention Time']

        print('Running alignment on ...')
                
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

                masterresults[c+'_se'][key] = std(masterresults[c][key]) / np.sqrt(masterresults['N Samples'][key])
                masterresults[c][key] = mean(masterresults[c][key])

        results_df = pd.DataFrame(masterresults).fillna(0)
        cols_at_end = [c for c in results_df.columns if 'Intensity' in c ]
        
        final_col_list = shared_columns + [ f for f in averaged_cols] + [ f + '_se' for f in averaged_cols] 

        final_col_list = [f for f in final_col_list if (f != 'file') & (f != 'Peak Height')] + used_elements + ['N Samples'] + cols_at_end
        
        results_df = results_df[final_col_list]
        
        return(results_df)
