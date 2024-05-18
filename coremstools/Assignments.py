__author__ = "Christian Dewey & Rene Boiteau"
__date__ = "2024 May 14"
__version__ = "0.0.1"

from pandas import DataFrame, read_csv

from coremstools.AssignmentError import AssignmentError
from coremstools.MolClassRetention import MolClassRetention
import coremstools.QualityControl as QualityControl
import coremstools.Dispersity as Dispersity
from coremstools.Parameters import Settings
import coremstools.Helpers as lcmsfns


class Assignments():
    """
    Base class for initial processing of CoreMS assignments. To be used before assembly of feature list. 

    Parameters 
    ----------
    sample_df : DataFrame 
        Pandas DataFrame containing a 'File' column with the name of each .raw file in the dataset (not full path). Defaults to None.
    t_interval : float
        Interval (min) over which scans are averaged in CoreMS LC assignments. Defaults to 2 (min).   
    """

    def __init__(self, sample_df, t_interval = 2):
        
        self.t_int = t_interval
        self.sample_df = sample_df     # dataframe

    def add_mol_class(self):
        '''
        Method to adds molecular class to CoreMS assignments & creates a 'Molecular Class' column. Altered assignment .csvs are rewritten with new columns. 
        '''
        for f in self.sample_df['File']:

            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            df = read_csv(fpath)
            heter = lcmsfns.get_heteroatoms(df)
            molclasses = lcmsfns.get_mol_class(heter)
            df2 = lcmsfns.assign_mol_class(df,molclasses)
            df2.to_csv(fpath, index = False)

    def run_internal_std_qc(self,timerange = [10,12]):
        '''
        Method to run the quality control checks with the internal standard m/z. Re-writes the sample list with additional columns for internal standard area, retention time, and QC pass/fail flag.

        Parameters
        ----------
        timerange : list
            A two-element list containing the minimum retention time and the maximum retention time between which to extract and plot the EIC for the internal standard. Defaults to [10, 12].
        
            '''

        self.sample_df = QualityControl.StandardQC(self.sample_df, timerange) ######################################### need to check why this returns something

    def run_assignment_error_plot(self, n_molclass = -1):
        '''
        Method to generate assignment error plots for QC checks. For each sample in the sample list, this method creates .jpg plots of (i) m/z Error (ppm) v. m/z and (ii) Molecular Classes of assignments over separation. The .jpgs are saved in the directory defined by Settings.assignments_directory.
        
        Parameters
        ----------
        n_molclass : int
            Specifies number of molecular classes to explicitly represent in error plots. If set to -1, all molecular classes will be explicitly represented. If set to a value greater than 0, the first n_molclass molecular classes, sorted from most abundant to least abundant, will be explicitly represented.
        '''

        print('plotting m/z error plots for ...')   
        for f in self.sample_df['File']:
            print('  '+ f)
            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            save_file = fpath.split('.')[0] + '_mz-error.jpg'
            AssignmentError.ErrorPlot(read_csv(fpath), save_file, n_molclass)
            
    
    def run_molclass_retention_plot(self, n_molclass = -1):
        '''
        Method to generate bar plots showing the proportions of molecular classes assigned to m/z within each time interval of the separation. Unassigned m/z are also shown.
        
        Parameters
        ----------
        n_molclass : int
            Specifies number of molecular classes to explicitly represent in plots. If set to -1, all molecular classes will be explicitly represented. If set to a value greater than 0, the first n_molclass molecular classes, sorted from most abundant to least abundant, will be explicitly represented. m/z with other molecular classes will be represented as 'Other'. Unassigned m/z are always represented. 
        '''
        print('plotting molecular classes v. retention time for ...')   
        for f in self.sample_df['File']:
            print('  '+ f)
            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            save_file = fpath.split('.')[0] + '_rt-mc.jpg'
            MolClassRetention.RTAssignPlot(read_csv(fpath), save_file, n_molclass)


    def run_dispersity_calculation(self):
        '''
        Method to runs dispersity calculation on each m/z in the CoreMS assignment file corresponding to each sample. The CoreMS assignment files are copied and saved as [SAMPLE_NAME] + dispersity_addend +'.csv' in the directory defined by Settings.assignments_directory. Currently quite slow. Would be good to do this calculation after the feature list is assembled.
        '''
        

        print('running dispersity calculation on ...')

        for f in self.sample_df['File']:
            print('  ' + f)
            fcsv = f.split('.')[0] + Settings.csvfile_addend + '.csv'
            Dispersity.CalculateDispersity(Settings.assignments_directory +  fcsv)






