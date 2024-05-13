from pandas import DataFrame

from coremstools.features.calc.Align import Align
from coremstools.features.calc.GapFill import GapFill 
from coremstools.Parameters import Settings

class Features(Align, GapFill):

    def __init__(self, sample_df):
        
        self.feature_list: DataFrame = None
        self.sample_df = sample_df

    def run_alignment(self):

        self.feature_list = Align.Align(self.sample_df)

    def run_gapfill(self):

        if self.feature_list is not None:
            self.feature_list = GapFill.GapFill(self.feature_list)
        else:
            self.run_alignment()
            self.feature_list = GapFill.GapFill(self.feature_list)

    def mz_error_flag(self):
        """
        This function identifies features with potentially significant mass measurement errors based on rolling average and standard deviation.

        Args:
            feature_list: A pandas DataFrame containing the aligned features with 'm/z Error (ppm)' and 'm/z Error (ppm) stdev' columns.

        Returns:
            The modified DataFrame with a new column 'mz error flag' indicating potential errors.
        """

        self.feature_list = self.feature_list.sort_values(by=['Calculated m/z'])
        self.feature_list['rolling error'] = self.feature_list['m/z Error (ppm)'].rolling(int(len(self.feature_list)/50), center=True,min_periods=0).mean()
        self.feature_list['mz error flag'] = abs(self.feature_list['rolling error'] - self.feature_list['m/z Error (ppm)']) / (4*self.feature_list['m/z Error (ppm) stdev'])

    def export(self):

        self.feature_list.to_csv(Settings.assignments_directory + 'feature_list.csv', index=False)

        
    def blank_flag(self):
        """
        This function calculates a 'blank' flag based on the intensity of a specific blank file compared to the maximum intensity in each feature's spectrum.

        Args:
            feature_list: A pandas DataFrame containing the aligned features with intensity columns prefixed with 'Intensity:'.
            blankfile: The filename of the blank data file (assumed to be defined elsewhere).

        Returns:
            The modified DataFrame with a new column 'blank' flag indicating potential blank contamination.
        """
        self.feature_list['Max Intense'] = self.feature_list.filter(regex='Intensity').max(axis=1)
        self.feature_list['blank'] = self.feature_list['Intensity:'+ blankfile.replace(dfiletype,'')].fillna(0)/feature_list['Max Intense']
        return(feature_list)