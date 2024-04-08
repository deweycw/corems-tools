from pandas import DataFrame, read_csv, unique
from numpy import mean, std, array, zeros, shape, where, log10
from tqdm import tqdm

class Features:


    def __init__(self, csv_list, include_dispersity = True):
        
        self.results: DataFrame = None
        self.shared_columns: list = None
        self.averaged_cols: list = None

        self.include_dispersity = include_dispersity
        self.csv_list = csv_list

        self._results = {}

        self._build_results_dict()



    def _build_results_dict(self):
        
        self.shared_columns = ['Time', 'Molecular Formula', 'Ion Charge', 'Calculated m/z', 'Heteroatom Class',  'DBE']

        self.averaged_cols = ['m/z',
                    'm/z Error (ppm)',
                    'Calibrated m/z',
                    'Resolving Power',
                    'm/z Error Score',
                    'Isotopologue Similarity',
                    'Confidence Score',
                    'S/N']
        
        if self.include_dispersity:

            self.averaged_cols.append('Dispersity')
            
        
        for c in self.shared_columns:
            
            self._results[c] = {}
        
        self._results['N Samples'] = {}
        
        for c in self.averaged_cols:
            
            self._results[c] = {}
            
            self._results[c + ' stdev'] = {}



    def Align(self):
        
        self.elements = []
        
        for file in self.csv_list:
            
            print('aligning ' + file)   

            result = read_csv(file)
            
            result = result[result['Molecular Formula'].notnull()]
            
            result['feature'] = list(zip(result['Time'],result['Molecular Formula']))
            
            file_name = file.replace('.csv','')
            
            self._results['Intensity: '+ file_name] = {}

            pbar = tqdm(range(len(result)), ncols=100)

            for ix in pbar:

                row = result.iloc[ix,:]

                if row['feature'] not in self._results['Time'].keys():

                    for c in self.shared_columns:

                        self._results[c][row['feature']] = row[c]
                    
                    current_elements = [x.rstrip('0123456789') for x in row['Molecular Formula'].split()]

                    for element in current_elements:

                        if element not in self.elements:

                            self.elements.append(element)

                            self._results[element] = {}

                        self._results[element][row['feature']] = row[element]
                    
                    self._results['Intensity: ' + file_name][row['feature']] = int(row['Peak Height'])

                    for c in self.averaged_cols:

                        self._results[c][row['feature']] = [row[c]]
                    
                else:
                    
                    self._results['Intensity: ' + file_name][row['feature']] = int(row['Peak Height'])

                    for c in self.averaged_cols:

                        self._results[c][row['feature']].append(row[c])


        print('writing N Samples column')

        pbar = tqdm(self._results['m/z'].keys(), ncols=100)

        for key in pbar:

            self._results['N Samples'][key] = len(self._results['m/z'][key])
            
            for c in self.averaged_cols:
                
                self._results[c + ' stdev'][key] = std(self._results[c][key])

                self._results[c][key] = mean(self._results[c][key])

            
        self.results = DataFrame(self._results).fillna(0)



    def LegacyGapFill(self):
        
        self.results['Gapfilled'] = False

        self.results['Low Conf. Gap'] = False 

        pbar = tqdm(range(len(self.results)))

        for i in pbar:

            row = self.results.iloc[i,:]

            resolution = row['Resolving Power']

            mass = row['Calibrated m/z']

            time = row['Time']

            mass_range = [mass* (1 - 2/resolution), mass * (1 + 2/resolution)]

            matches = self.results[(self.results['Calibrated m/z'] > mass_range[0]) &
                                   (self.results['Calibrated m/z'] < mass_range[1]) &
                                   (self.results['Time'] == time)]
            
            if (len(matches) > 1):

                self.results.loc[i,'Gapfilled'] = True

                results_intensity_cols = self.results.filter(regex='Intensity').columns

                self.results.loc[i,results_intensity_cols] = matches.filter(regex='Intensity').sum(axis=0)
                
                #print(self.results.loc[i,results_intensity_cols])

                if row['Confidence Score'] < max(matches['Confidence Score']):

                    self.results.loc[i,'Low Conf. Gap'] = True



    def GapFill(self):

        print('performing gap fill')

        intensity_cols = list(self.results.filter(regex='Intensity').columns)

        cols = self.results.columns

        self.results = self.results[[col for col in cols if 'Intensity' not in col] + [col for col in cols if 'Intensity' in col]]

        self.results.sort_values(['Time','Calibrated m/z'], inplace=True)

        self.results['Multiple Peaks w/in Uncertainty'] = None

        self.results['Gapfill ID'] = None

        self.results['Gapfill Flag'] = None

        self.results['Gapfill Molecular Formula'] = None

        for time_step in unique(self.results['Time']):

            time_step_df = self.results[self.results['Time'] == time_step].copy()

            n_rows = len(time_step_df)

            mz_array = array([time_step_df['Calibrated m/z'] for i in range(n_rows)])

            R_array = array([time_step_df['Resolving Power'] for i in range(n_rows)])

            FWHM_array = mz_array / R_array

            mz_error_array = FWHM_array + FWHM_array.T

            mz_diff_array = abs(mz_array - mz_array.T)

            gapfill_inds_array = array(mz_diff_array < mz_error_array)

            n_gapfills_inds = where(gapfill_inds_array)

            n_gapfills_array = zeros(shape(mz_array))

            n_gapfills_array[n_gapfills_inds[0], n_gapfills_inds[1]] = 1

            gapfill_sum = n_gapfills_array.sum(axis=0)

            n_gapfills_vector = array([True if i > 1 else None for i in gapfill_sum ])
            
            time_step_df['Multiple Peaks w/in Uncertainty'] = n_gapfills_vector
            
            time_step_df.sort_values(['Time','Calibrated m/z'], inplace=True)

            offset_diag_rows = array([ i for i in range(0,shape(mz_array)[0]-1)])

            offset_diag_cols = array([ j for j in range(1, shape(mz_array)[0])])

            neighboring_mz_diffs = mz_diff_array[offset_diag_rows, offset_diag_cols]
            
            neighboring_mz_err = mz_error_array[offset_diag_rows, offset_diag_cols]

            residual_diff = neighboring_mz_diffs - neighboring_mz_err

            transition_inds = where(residual_diff > 0 )

            n_true_block = 1
                        
            gapfill_column = zeros((shape(mz_array)[0],1))

            pbar = tqdm(range(len(time_step_df)),ncols=100)

            for ix in pbar:
                
                pbar.set_description_str(desc="Adding gapfill ID for timestep %s" %(time_step) , refresh=True)

                gapfill_bool = time_step_df.iloc[ix].loc['Multiple Peaks w/in Uncertainty']

                if gapfill_bool == True:
                    
                    if log10(n_true_block) < 1:

                        add_string = '.00'
                    
                    elif (log10(n_true_block) >= 1) &(log10(n_true_block) < 2):

                        add_string = '.0'

                    else:

                        add_string = '.'

                    gap_id = float(str(time_step) + add_string + str(n_true_block))
                    
                    gapfill_column[ix,0] = gap_id

                    if ix in transition_inds[0]:

                        n_true_block = n_true_block + 1 

            time_step_df['Gapfill ID'] = gapfill_column
            
            gap_ids_list = [id for id in unique(time_step_df['Gapfill ID']) if id > 0]


            pbar = tqdm(gap_ids_list, ncols=100)
            
            for id in pbar:

                pbar.set_description_str(desc="Adding gapfill flag for timestep %s" %(time_step) , refresh=True)

                id_df = time_step_df[time_step_df['Gapfill ID'] == id].copy()

                id_intensity_sum = id_df.filter(regex='Intensity').sum(axis=0)

                id_df[intensity_cols] = id_intensity_sum

                id_max_confidence_score = max(id_df['Confidence Score'])

                id_df.loc[id_df['Confidence Score'] < id_max_confidence_score,'Gapfill Flag'] = True

                gapfilled_mf = id_df.loc[id_df['Confidence Score'] == id_max_confidence_score,'Molecular Formula']

                id_df.loc[id_df['Confidence Score'] < id_max_confidence_score,'Gapfill Molecular Formula'] = gapfilled_mf.iloc[0]

                time_step_df[time_step_df['Gapfill ID'] == id] = id_df

            
            
            
            self.results[self.results['Time'] == time_step] = time_step_df

            
        