from pandas import DataFrame, read_csv
from numpy import mean, std
from tqdm import tqdm

from coremstools.Parameters import Settings

def Align(sample_df):

    def build_masterresults_dict(shared_columns, averaged_cols):
        
        masterresults={}

        for col in shared_columns:
            
            masterresults[col] = {}
        
        masterresults['N Samples'] = {}
        
        for col in averaged_cols:
            
            masterresults[col] = {}
            
            masterresults[col + ' stdev'] = {}

        return masterresults

    assignments_dir = Settings.assignments_directory
    disp_addend = Settings.dispersity_addend

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
    print('running alignment on ...')
            
    elements=[]

    masterresults = build_masterresults_dict(shared_columns, averaged_cols)

    for file in sample_df['File']:

        print('  ' + file)

        file = assignments_dir + file.split('.')[0] + disp_addend + '.csv'

        results = read_csv(file)
        
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
    

