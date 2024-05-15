__author__ = "Christian Dewey & Rene Boiteau"
__date__ = "2024 May 14"
__version__ = "0.0.1"

from pandas import DataFrame, read_csv

import coremstools.AssignmentError as AssignmentError
import coremstools.QualityControl as QualityControl
import coremstools.Dispersity as Dispersity
from coremstools.Parameters import Settings
import coremstools.Helpers as lcmsfns


class Assignments():
    """
    Base class for initial processing of CoreMS assignments. Used before assembly of feature list. 

    Parameters 
    ----------
    sample_df : DataFrame 
        Pandas DataFrame containing a 'File' column with the name of each .raw file in the dataset (not full path). Defaults to None.
    t_interval : float
        Interval (min) over which scans are averaged in CoreMS LC assignments. Defaults to 2 (min).   

    Methods
    -------
    add_mol_class()
        Adds molecular class to CoreMS assignments & creates a 'Molecular Class' column. Altered assignment .csvs are rewritten with new columns. 
    run_internal_std_qc(timerange=[10,12]) -> DataFrame
        Runs the quality control checks with internal standard m/z. Returns copy of sample list DataFrame with additional columns for internal standard area, retention time, and QC pass/fail flag.
    run_assignment_error_plot()
        For each sample in the sample list, this method creates .jpg plots of (i) m/z Error (ppm) v. m/z and (ii) Molecular Classes of assignments over separation. The .jpgs are saved in the directory defined by Settings.assignments_directory.
    run_dispersity_calculation()
        Runs dispersity calculation on each m/z in the CoreMS assignment file corresponding to each sample. The CoreMS assignment files are copied and saved as [SAMPLE_NAME] + dispersity_addend +'.csv' in the directory defined by Settings.assignments_directory. ***** Currently quite slow. Would be good to do this calculation after the feature list is assembled. ******
    """

    def __init__(self, sample_df, t_interval = 2):
        
        self.t_int = t_interval
        self.sample_df = sample_df     # dataframe

    def add_mol_class(self):
        for f in self.sample_df['File']:

            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            df = read_csv(fpath)
            heter = lcmsfns.get_heteroatoms(df)
            molclasses = lcmsfns.get_mol_class(heter)
            df2 = lcmsfns.assign_mol_class(df,molclasses)
            df2.to_csv(fpath, index = False)

    def run_internal_std_qc(self,timerange=[10,12]):

        self.sample_df = QualityControl.StandardQC(self.sample_df, timerange) ######################################### need to check why this returns something

    def run_assignment_error_plot(self):

        for f in self.sample_df['File']:

            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            save_file = fpath.split('.')[0] + '_mz-error.jpg'
            AssignmentError.ErrorPlot(read_csv(fpath), save_file)
            save_file = fpath.split('.')[0] + '_rt-error.jpg'
            AssignmentError.RTAssignPlot(read_csv(fpath), save_file)

    def run_dispersity_calculation(self):

        print('running dispersity calculation on ...')

        for f in self.sample_df['File']:
            print('  ' + f)
            Dispersity.CalculateDispersity(f, self.t_int)






