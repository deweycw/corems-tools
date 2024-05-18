from pandas import DataFrame
import coremstools.Align as Align
import coremstools.GapFill as GapFill 
from coremstools.Parameters import Settings
import coremstools.Dispersity as Dispersity


class Features:
    """
    Base class for holding CoreMS features across a dataset. 

    Parameters 
    ----------
    sample_list : DataFrame 
        Pandas DataFrame containing a 'File' column with the name of each .raw file in the dataset (not full path). Defaults to None.

    Methods
    -------
    run_alignment()
        Aligns features across dataset. 
    run_gapfill()
        Runs gapfill. 
    flag_errors()
        Identifies features with potentially significant mass measurement errors based on rolling average and standard deviation of measured m/z.
    flag_blank_features()
        Calculates a 'blank' flag based on the intensity of a specific blank file compared to the maximum intensity in each feature's spectrum.
    export()
        Writes feature list to .csv file. 
        
    """
    def __init__(self, sample_list):
        
        self.feature_list: DataFrame = None
        self.sample_list = sample_list

    def run_alignment(self):

        self.feature_list = Align.Align(self.sample_list)

    def run_gapfill(self):

        if self.feature_list is not None:
            self.feature_list = GapFill.GapFill(self.feature_list)
        else:
            self.run_alignment()
            self.feature_list = GapFill.GapFill(self.feature_list)

    def flag_errors(self):

        self.feature_list = self.feature_list.sort_values(by=['Calculated m/z'])
        self.feature_list['rolling error'] = self.feature_list['m/z Error (ppm)'].rolling(int(len(self.feature_list)/50), center=True,min_periods=0).mean()
        self.feature_list['mz error flag'] = abs(self.feature_list['rolling error'] - self.feature_list['m/z Error (ppm)']) / (4*self.feature_list['m/z Error (ppm) stdev'])

    def export(self):

        self.feature_list.to_csv(Settings.assignments_directory + 'feature_list.csv', index=False)

        
    def flag_blank_features(self, blank_sample):

        col = None

        for col in self.feature_list.columns:
            
            if blank_sample.split('.')[0] in col:

                blank_sample_col = col

        self.feature_list['Max Intensity'] = self.feature_list.filter(regex='Intensity').max(axis=1)
        self.feature_list['blank'] = self.feature_list[blank_sample_col].fillna(0) / self.feature_list['Max Intensity']
    


    def run_dispersity_calculation(self):
        '''
        Method to run dispersity calculation on each feature list. The CoreMS assignment files are copied and saved as [SAMPLE_NAME] + dispersity_addend +'.csv' in the directory defined by Settings.assignments_directory. 
        '''
        
        print('running dispersity calculation on ...')

        for f in self.sample_list['File']:
            print('  ' + f)
            Dispersity.CalculateDispersity(f)

class _FeatureListDispersity(Features):

    def __init__(self, feature_list,sample_list):
        super.__init__(feature_list, sample_list)