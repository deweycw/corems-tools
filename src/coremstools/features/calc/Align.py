from pandas import DataFrame, read_csv
from numpy import mean, std
from tqdm import tqdm

class Align:

    def Align(self):
        
        self.elements = []
        
        for file in self.csv_list:
            
            print('aligning ' + file)   

            result = read_csv(file)
            
            result = result[result['Molecular Formula'].notnull()]
            
            result['feature'] = list(zip(result['Time'],result['Molecular Formula']))
            
            file_name = file.replace('.csv','').split('/')[-1]
            
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


