import warnings
warnings.filterwarnings('ignore')

from coremstools.Parameters import Settings
from coremstools.DataSet import DataSet

if __name__ == '__main__':
    
    Settings.raw_file_directory = '/home/christiandewey/Documents/corems/rawfiles/done/'
    Settings.assignments_directory = Settings.raw_file_directory
    Settings.internal_std_mz = 678.2915
    Settings.std_time_range = [10,12]
    Settings.time_interval = 1
    #Settings.blank_sample_name = '20221103_LBA_Boiteau_Zorbax3p5_qh2o_fullmz'
    Settings.csvfile_addend = '_cd-2405'
    Settings.dispersity_addend = '_dispersity'

    dset = DataSet(path_to_sample_list = Settings.raw_file_directory + 'SampleList.csv')
    
    #dset.run_internal_std_qc()

    #dset.run_assignment_error_plots()

    #dset.run_molclass_retention_plots()

    #dset.run_dispersity_calcs()

    dset.run_alignment()

    #dset.export_feature_list()

    #dset.run_gapfill()

    #dset.flag_blank_features()

    

    

