from pandas import DataFrame

from coremstools.features.calc.Align import Align
from coremstools.features.calc.GapFill import GapFill 

class Features(Align, GapFill):

    def __init__(self, csv_list, include_dispersity = True):
        
        self.results: DataFrame = None
        self.shared_columns: list = None
        self.averaged_cols: list = None

        self.include_dispersity = include_dispersity
        self.csv_list = csv_list

        self._results = {}

        self._build_results_dict()


    def _build_results_dict(self):
        
        self.shared_columns = ['Time', 'Molecular Formula','mol_class', 'Ion Charge', 'Calculated m/z', 'Heteroatom Class',  'DBE']

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

