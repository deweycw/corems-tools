from pandas import DataFrame, read_csv
from numpy import mean, std
from tqdm import tqdm

class Features:


    def __init__(self):
        
        self.results: DataFrame = None
        self.shared_columns: list = None
        self.averaged_cols: list = None

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
                    #'Dispersity']
        
        
        for c in self.shared_columns:
            
            self._results[c] = {}
        
        self._results['N Samples'] = {}
        
        for c in self.averaged_cols:
            
            self._results[c] = {}
            
            self._results[c + ' stdev'] = {}



    def Align(self, csv_list):
        
        self.elements = []
        
        for file in csv_list:
            
            print('aligning ' + file)   

            result = read_csv(file)
            
            result = result[result['Molecular Formula'].notnull()]
            
            result['feature'] = list(zip(result['Time'],result['Molecular Formula']))
            
            file_name = file.replace('.csv','')
            
            self._results['Intensity: '+ file_name] = {}

            pbar = tqdm(range(len(result)))

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

        pbar = tqdm(self._results['m/z'].keys())

        for key in pbar:

            self._results['N Samples'][key] = len(self._results['m/z'][key])
            
            for c in self.averaged_cols:
                
                self._results[c + ' stdev'][key] = std(self._results[c][key])

                self._results[c][key] = mean(self._results[c][key])

            
        self.results = DataFrame(self._results).fillna(0)



    def GapFill(self):
        
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

    def VGapFill(self):

        from numpy import array, zeros, shape, where, any
        from pandas import unique

        self.results['>1 Peak w/in Uncertainty'] = None

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

            n_gapfills_vector = array([True if i > 1 else False for i in gapfill_sum ])
            
            time_step_df['>1 Peak w/in Uncertainty'] = n_gapfills_vector
            
            time_step_df.sort_values(['Time','Calibrated m/z'], inplace=True)

            self.results.sort_values(['Time','Calibrated m/z'], inplace=True)
            
            self.results[self.results['Time'] == time_step] = time_step_df





