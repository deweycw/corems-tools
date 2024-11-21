__author__ = "Christian Dewey & Rene Boiteau"
__date__ = "2024 May 14"
__version__ = "0.0.1"

from pandas import read_csv

from coremstools.FeatureList import Features
from coremstools.Parameters import Settings
from coremstools.QualityControl import QualityControl
from coremstools.AssignmentError import AssignmentError
from coremstools.MolClassRetention import MolClassRetention
from coremstools.Dispersity import Dispersity

import coremstools.AssignmentCalcs as lcmsfns


class DataSet(Features):
    """
    Base class for CoreMS dataset object.

    Parameters 
    ----------
    path_to_sample_list : str 
        Full path to dataset .csv file. Must contain a column called 'File' with a list of the Thermo .raw files in the dataset. Required if sample_list is not defined. 
    sample_list : DataFrame
        DataFrame containing sample list. Must contain a column called 'File' with a list of the Thermo .raw files in the dataset. Required if path_to_sample_list is not defined. 

    """

    def __init__(self, path_to_sample_list=None, sample_list=None):
        
        self.path_to_sample_list = path_to_sample_list  
        self.sample_list = sample_list

        self.time_interval = Settings.time_interval
        self.feature_list = None
        self.feature_list_df = None

        if (self.sample_list == None) & (self.path_to_sample_list != None):
            
            self._initialize_from_sample_list_file()
        
        elif (self.sample_list == None) & (self.path_to_sample_list == None):

            print('Please provide either (1) the sample list as a Pandas DataFrame or (2) a path to sample list (.csv)')

        elif (self.sample_list != None) & (self.path_to_sample_list != None):

            print('A sample list dataframe and a path were provived. Defaulting to provided DataFrame.')
    

    def _initialize_from_sample_list_file(self):

        self.sample_list = read_csv(self.path_to_sample_list, index_col = None)

    def assign_mol_class(self):
        '''
        Method to assign Molecular Classes to CoreMS assignment results.
        ''' 

        for f in self.sample_list['File']:
            df = read_csv(Settings.assignments_directory + f.replace('.raw','.csv'))
            df = lcmsfns.add_mol_class(df)
            df.to_csv(Settings.assignments_directory + f.replace('.raw','.csv'))


    def run_internal_std_qc(self):
        '''
        Method to run the quality control checks with the internal standard m/z for all samples in dataset. Re-writes the sample list with additional columns for internal standard area, retention time, and QC pass/fail flag.
        '''

        self.sample_list = QualityControl.StandardQC(self, self.sample_list)

        if self.path_to_sample_list == None:
            self.path_to_sample_list == Settings.assignments_directory + 'sample_list.csv'
        
        print('Sample list saved to assignments directory, with pass/fail columns')
        self.sample_list.to_csv(self.path_to_sample_list, index =False)

    def _check_for_molclass(self, assignments, fpath):

        if 'Molecular Class' not in assignments.columns:
            assignments = lcmsfns.add_mol_class(assignments)
            assignments.to_csv(fpath, index = False)


    def run_assignment_error_plots(self, n_molclass = -1):
        '''
        Method to generate assignment error plots for QC checks. For each sample in the sample list, this method creates .jpg plots of (i) m/z Error (ppm) v. m/z and (ii) Molecular Classes of assignments over separation. The .jpgs are saved in the directory defined by Settings.assignments_directory.
        
        Parameters
        ----------
        n_molclass : int
            Specifies number of molecular classes to explicitly represent in error plots. If set to -1, all molecular classes will be explicitly represented. If set to a value greater than 0, the first n_molclass molecular classes, sorted from most abundant to least abundant, will be explicitly represented.
        '''

        print('\nplotting m/z error plots for ...')   
        for f in self.sample_list['File']:
            print('  '+ f)
            fpath = Settings.assignments_directory + f.split('.')[0] + '.csv'
            save_file = fpath.split('.')[0] + '_mz-error.jpg'

            assignments = read_csv(fpath)

            self._check_for_molclass(assignments, fpath)

            AssignmentError.ErrorPlot(self, assignments, save_file, n_molclass)
            
    
    def run_molclass_retention_plots(self, n_molclass = -1):
        '''
        Method to generate bar plots showing the proportions of molecular classes assigned to m/z within each time interval of the separation. Unassigned m/z are also shown.
        
        Parameters
        ----------
        n_molclass : int
            Specifies number of molecular classes to explicitly represent in plots. If set to -1, all molecular classes will be explicitly represented. If set to a value greater than 0, the first n_molclass molecular classes, sorted from most abundant to least abundant, will be explicitly represented. m/z with other molecular classes will be represented as 'Other'. Unassigned m/z are always represented. 
        '''
        print('\nplotting molecular classes v. retention time for ...')   
        for f in self.sample_list['File']:
            print('  '+ f)
            fpath = Settings.assignments_directory + f.split('.')[0]+ '.csv'

            assignments = read_csv(fpath)

            self._check_for_molclass(assignments, fpath)

            save_file = fpath.split('.')[0] + '_rt-mc.jpg'
            MolClassRetention.RTAssignPlot(self, assignments, save_file, n_molclass)


    def run_dispersity_calcs(self):
        '''
        Method to runs dispersity calculation on each m/z in the CoreMS assignment file corresponding to each sample. The CoreMS assignment files are copied and saved as [SAMPLE_NAME] + '.csv' in the directory defined by Settings.assignments_directory. Currently quite slow. Would be good to do this calculation after the feature list is assembled.
        '''

        print('\nrunning dispersity calculation on ...')

        for f in self.sample_list['File']:
            print('  ' + f)
            fcsv = f.split('.')[0] + '.csv'
            Dispersity.CalculateDispersity(Settings.assignments_directory +  fcsv)

    def _check_for_feature_list(self):

        if self.feature_list == None:
            self.feature_list = Features(self.sample_list)
            try:
                if len(self.feature_list_df.columns) > 0:
                    self.feature_list.feature_list_df = self.feature_list_df
            except:
                pass
        

    def run_alignment(self, include_dispersity = True, experimental = False):
        """
        Method to assemble an aligned feature list for the dataset. The aligned feature list is a dataframe containing a row for each [molecular formula]-[retention time] pair (what we call a feature) in the entire dataset. The dataframe contains the intensity of each feature in each sample in the data, as well as the average and stdev of each of the following parameters: measured m/z of the feature; calibrated m/z of the feature; resolving power of the instrument at the measured m/z; m/z error score; istopologue similarity score; confidence score; S/N; and dispersity. 

        Parameters
        ----------
        include_dispersity : bool
            Flag indicating whether or not to include dispersity in feature list. If set to False, dispersity will not be include. Default is True. 
        
        experimental : bool
            Flag indicating whether to use experimental funcationality for alignment. Default is False.
        """

        self._check_for_feature_list()
        self.feature_list.run_alignment(include_dispersity, experimental)


    def run_consolidation(self, gapfill_variable = 'Confidence Score', include_dispersity = True, experimental = False):
        '''
        Method to perform gapfilling of features across dataset.
        '''        
        self._check_for_feature_list()
        self.feature_list.run_consolidation(gapfill_variable, include_dispersity, experimental)


    def calc_rolling_error(self):

        self._check_for_feature_list()
        self.feature_list.flag_errors()


    def flag_blank_features(self):
        '''
        Method to flag features that appear in blank.
        '''
        self._check_for_feature_list()
        self.feature_list.flag_blank_features()


    def export_feature_list(self, fname = 'feature_list.csv'):
        '''
        Method to export feature list as .csv file. Will be exported to the directory defined in Parameters.Settings.assignments_directory.
        '''
        self._check_for_feature_list()
        self.feature_list.export_csv(fname)
        


