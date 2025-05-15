#from pandas import DataFrame
#import dask.dataframe as dd

from coremstools.Align import Align
from coremstools.Consolidate import Consolidate 
from coremstools.Parameters import Settings
import numpy as np

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
        Identifies features with potentially significant mass measurement errors based on rolling average and standard error of measured m/z.
    flag_blank_features()
        Calculates a 'blank' flag based on the intensity of a specific blank file compared to the maximum intensity in each feature's spectrum.
    export()
        Writes feature list to .csv file. 
    """
    def __init__(self, sample_list):
        
        self.feature_list_df = None
        self.sample_list = sample_list

    def run_alignment(self, include_dispersity, experimental):
        
        if experimental:
        
            self.feature_list_df = Align.Align_exp(self, self.sample_list, include_dispersity)

        else:

            self.feature_list_df = Align.run(self, self.sample_list, include_dispersity)


    def run_consolidation(self, gapfill_variable, include_dispersity):

        if self.feature_list_df is not None:         
            self.feature_list_df = Consolidate.run(self, gapfill_variable, self.feature_list_df)

        else:
            self.run_alignment(include_dispersity)
            self.feature_list_df = Consolidate.run(self, self.feature_list_df)
        

    def flag_errors(self, n_iter=3):

        '''
        Method that (1) calculates a rolling average of the assignment error, from lowest to highest calculated m/z, for each feature in the feature list, and (2) calculates an error flag, which is the absolute value of the difference of the rolling average error and the average error of the individual feature divided by 4 times the standard deviation of the m/z error for the feature. 
        '''
        print('Running holistic m/z error filter...')

        self.feature_list_df = self.feature_list_df.sort_values(by=['Calculated m/z'])
        if 'consolidated flag' in self.feature_list_df.columns:
            self.feature_list_df['rolling error'] = self.feature_list_df[self.feature_list_df['consolidated flag'] == 0]['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()
            self.feature_list_df['rolling error'].interpolate(method='linear', inplace=True)

        else:
            self.feature_list_df['rolling error'] = self.feature_list_df['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()

        self.feature_list_df['mz error flag'] = abs(self.feature_list_df['rolling error'] - self.feature_list_df['m/z Error (ppm)']) / (4*self.feature_list_df['m/z Error (ppm)_se'])

        if n_iter>1:
            for i in range(n_iter-1):
                self.feature_list_df['rolling error'] = np.nan
                if 'consolidated flag' in self.feature_list_df.columns:
                    self.feature_list_df['rolling error'] = self.feature_list_df[(self.feature_list_df['consolidated flag'] == 0) & (self.feature_list_df['mz error flag']<1)]['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()
                    self.feature_list_df['rolling error'].interpolate(method='linear', inplace=True)

                else:
                    self.feature_list_df['rolling error'] = self.feature_list_df[self.feature_list_df['mz error flag']<1]['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()

                #self.feature_list_df['mz error flag'] = abs(self.feature_list_df['rolling error'] - self.feature_list_df['m/z Error (ppm)']) / (4*self.feature_list_df['m/z Error (ppm)_se']*np.sqrt(self.feature_list_df['N Samples']))
                self.feature_list_df['mz error flag'] = abs(self.feature_list_df['rolling error'] - self.feature_list_df['m/z Error (ppm)']) / (4*self.feature_list_df['m/z Error (ppm)_se'])


    def flag_blank_features(self):
        """
        This function calculates and adds to a feature the following columns:
            'Max Intensity': The maximum signal intensity observed across all samples
            'Max Blank': The maximum signal intensity observed across blank files.
            'Blank' column based on the max intensity of a specific set of blank file compared to the maximum intensity in each feature's spectrum.

        Args:
            self.feature_list_df: A pandas DataFrame containing the aligned features with intensity columns prefixed with 'Intensity:'.
            Settings.blank_sample_list: A list of the filename of the blank data file.


        """

        print('Flagging blank features...')

        if self.feature_list_df is None:

            self.run_alignment()

        blank_col = []

        for blank_sample in Settings.blank_sample_list:

            if '.' in blank_sample:
            
                blank_sample = blank_sample.split('.')[0]


            for col in self.feature_list_df.columns:
                
                if blank_sample in col:

                    blank_col.append(col)

        
        #blank_col = ["Intensity:" + item for item in Settings.blank_sample_list]
        self.feature_list_df['Max Intensity']=self.feature_list_df.filter(regex='Intensity').max(axis=1)
        self.feature_list_df['Max Blank']=self.feature_list_df[blank_col].max(axis=1)
        self.feature_list_df['blank']=self.feature_list_df['Max Blank']/self.feature_list_df['Max Intensity']



    def stoichiometric_classification(self):

        print('Determining stoichiometric classifications...')

        self.feature_list_df['Stoichiometric classification']='Unclassified'

        self.feature_list_df['O/C']=self.feature_list_df['O']/self.feature_list_df['C']
        self.feature_list_df['H/C']=self.feature_list_df['H']/self.feature_list_df['C']
        # Calculate atomic stoichiometries
        contains_N = True
        contains_P = True
        cols_to_remove = []
        if not 'N' in self.feature_list_df.columns:
            self.feature_list_df['N']=0
            cols_to_remove = cols_to_remove + ['N','N/C']
            contains_N = False
        if not 'P' in self.feature_list_df.columns:
            self.feature_list_df['P']=0
            cols_to_remove = cols_to_remove + ['P', 'P/C']
            contains_P = False
        if not 'S' in self.feature_list_df.columns:
            self.feature_list_df['S']=0
            cols_to_remove.append('S')

        if (not contains_N) or (not contains_P):
            cols_to_remove.append('N/P')
        

        self.feature_list_df['N/C']=self.feature_list_df['N']/self.feature_list_df['C']
        self.feature_list_df['P/C']=self.feature_list_df['P']/self.feature_list_df['C']
        self.feature_list_df['N/P']=self.feature_list_df['N']/self.feature_list_df['P']

        self.feature_list_df['NOSC'] =  4 -(4*self.feature_list_df['C'] 
                                + self.feature_list_df['H'] 
                                - 3*self.feature_list_df['N'] 
                                - 2*self.feature_list_df['O'])/self.feature_list_df['C']

        self.feature_list_df.loc[(self.feature_list_df['O/C']<=0.6) & 
                            (self.feature_list_df['H/C']>=1.32) & 
                            (self.feature_list_df['N/C']<=0.126) &
                            (self.feature_list_df['P/C']<0.35)
                            ,'Stoichiometric classification'] = 'Lipid'

        self.feature_list_df.loc[(self.feature_list_df['O/C']<=0.6) & 
                            (self.feature_list_df['H/C']>=1.32) & 
                            (self.feature_list_df['N/C']<=0.126) &
                            (self.feature_list_df['P/C']<0.35) &
                            (self.feature_list_df['P']>0)
                            ,'Stoichiometric classification'] = 'Phospholipid'

        self.feature_list_df.loc[(self.feature_list_df['O/C']>=0.61) & 
                            (self.feature_list_df['H/C']>=1.45) & 
                            (self.feature_list_df['N/C']>0.07) & 
                            (self.feature_list_df['N/C']<=0.2) & 
                            (self.feature_list_df['P/C']<0.3) & 
                            (self.feature_list_df['O']>=3) &
                            (self.feature_list_df['N']>=1)
                            ,'Stoichiometric classification'] = 'A-Sugar'

        self.feature_list_df.loc[(self.feature_list_df['O/C']>=0.8) & 
                            (self.feature_list_df['H/C']>=1.65) & 
                            (self.feature_list_df['H/C']<2.7) &
                            (self.feature_list_df['O']>=3) &
                            (self.feature_list_df['N']==0)
                            ,'Stoichiometric classification'] = 'Carbohydrate'

        self.feature_list_df.loc[(self.feature_list_df['O/C']>=0.5) & 
                            (self.feature_list_df['O/C']<1.7) & 
                            (self.feature_list_df['H/C']>1) & 
                            (self.feature_list_df['H/C']<1.8) &
                            (self.feature_list_df['N/C']>=0.2) & 
                            (self.feature_list_df['N/C']<=0.5) & 
                            (self.feature_list_df['N']>=2) &
                            (self.feature_list_df['P']>=1) &
                            (self.feature_list_df['S']==0) &
                            (self.feature_list_df['Calculated m/z']>305) &
                            (self.feature_list_df['Calculated m/z']<523)
                            ,'Stoichiometric classification'] = 'Nucleotide'

        self.feature_list_df.loc[(self.feature_list_df['O/C']<=1.15) & 
                            (self.feature_list_df['H/C']<1.32) & 
                            (self.feature_list_df['N/C']<0.126) &
                            (self.feature_list_df['P/C']<=0.2) 
                            ,'Stoichiometric classification'] = 'Phytochemical'

        self.feature_list_df.loc[(self.feature_list_df['S']>0)
                            ,'Stoichiometric classification'] = 'Organosulfur'

        self.feature_list_df.loc[(self.feature_list_df['O/C']>0.12) & 
                            (self.feature_list_df['O/C']<=0.6) & 
                            (self.feature_list_df['H/C']>0.9) & 
                            (self.feature_list_df['H/C']<2.5) & 
                            (self.feature_list_df['N/C']>=0.126) & 
                            (self.feature_list_df['N/C']<=0.7) & 
                            (self.feature_list_df['P/C']<0.17) & 
                            (self.feature_list_df['N']>=1)
                            ,'Stoichiometric classification'] = 'Peptide'

        self.feature_list_df.loc[(self.feature_list_df['O/C']>0.6) & 
                            (self.feature_list_df['O/C']<=1) & 
                            (self.feature_list_df['H/C']>1.2) & 
                            (self.feature_list_df['H/C']<2.5) & 
                            (self.feature_list_df['N/C']>=0.2) & 
                            (self.feature_list_df['N/C']<=0.7) & 
                            (self.feature_list_df['P/C']<0.17) & 
                            (self.feature_list_df['N']>=1)
                            ,'Stoichiometric classification'] = 'Peptide'

        self.feature_list_df.loc[(self.feature_list_df['Is Isotopologue']>0),'Stoichiometric classification']='Isotoplogue'
        for col in cols_to_remove:
            self.feature_list_df.drop(col, axis = 1, inplace=True)

    def export_csv(self, fname):

        print('writing to .csv...')
        dir = Settings.assignments_directory
        self.feature_list_df.to_csv(dir + fname, index = False) #, single_file = True, header_first_partition_only = True)
        

    def export_parquet(self):

        print('writing to .parquet...')
        dir = '/home/christiandewey/Dropbox/'
        #dir = Settings.assignments_directory
        self.feature_list_df.to_parquet(dir + 'feature_list.parquet', compute =True)
        