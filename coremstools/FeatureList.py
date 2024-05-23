from pandas import DataFrame
import dask.dataframe as dd

from coremstools.Align import Align
from coremstools.GapFill import GapFill 
from coremstools.Parameters import Settings

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
        
        self.feature_list_ddf = None
        self.sample_list = sample_list

    def run_alignment(self,experimental = False):
        
        if experimental:
        
            self.feature_list_ddf = Align.Align_exp(self, self.sample_list)

        else:

            self.feature_list_ddf = Align.Align(self, self.sample_list)


    def run_gapfill(self):

        if self.feature_list_ddf is not None:
            self.feature_list_ddf = GapFill.GapFill(self, self.feature_list_ddf)
        else:
            self.run_alignment()
            self.feature_list_ddf = GapFill.GapFill(self, self.feature_list_ddf)
        
    def flag_errors(self):

        '''
        Method that (1) calculates a rolling average of the assignment error, from lowest to highest calculated m/z, for each feature in the feature list, and (2) calculates an error flag, which is the absolute value of the difference of the rolling average error and the average error of the individual feature divided by 4 times the standard deviation of the m/z error for the feature. '''

        self.feature_list_df.sort_values(by=['Calculated m/z'], inplace=True)

        self.feature_list_df['rolling error'] = self.feature_list_df['m/z Error (ppm)'].rolling(int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()

        self.feature_list_df['mz error flag'] = abs(self.feature_list_df['rolling error'] - self.feature_list_df['m/z Error (ppm)']) / (4*self.feature_list_df['m/z Error (ppm) stdev'])

        
    def flag_blank_features(self):

        print('flagging blank features')

        if self.feature_list_df is None:

            self.run_alignment()

        col = None

        blank_sample = Settings.blank_sample_name

        if '.' in blank_sample:
        
            blank_sample = blank_sample.split('.')[0]

        for col in self.feature_list_df.columns:
            
            if blank_sample in col:

                blank_sample_col = col

        self.feature_list_df['Max Intensity'] = self.feature_list_df.filter(regex='Intensity').max(axis=1)
        self.feature_list_df['blank'] = self.feature_list_df[blank_sample_col].fillna(0) / self.feature_list_df['Max Intensity']

        #self.feature_list_df.to_csv(Settings.assignments_directory + 'feature_list_df.csv')


    def export_csv(self):

        print('writing to .csv...')
        dir = '/home/christiandewey/Dropbox/'
        #dir = Settings.assignments_directory
        self.feature_list_ddf.to_csv(dir + 'feature_list.csv',index = False) #, single_file = True, header_first_partition_only = True)
        
        #self.feature_list_ddf.to_csv(Settings.assignments_directory + 'feature_list-leg.csv',index = False)

    def export_parquet(self):

        print('writing to .parquet...')
        dir = '/home/christiandewey/Dropbox/'
        #dir = Settings.assignments_directory
        self.feature_list_ddf.to_parquet(dir + 'feature_list.parquet', compute =True)
        
        #self.feature_list_ddf.to_csv(Settings.assignments_directory + 'feature_list-leg.csv',index = False)