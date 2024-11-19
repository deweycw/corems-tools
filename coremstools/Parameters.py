__author__ = "Christian Dewey & Rene Boiteau"
__date__ = "2024 May 14"
__version__ = "0.0.1"

import dataclasses

@dataclasses.dataclass
class Settings:
    """
    Settings class to store global settings for processing CoreMS assignments.

    Attributes
    ----------
    raw_file_directory : str
        Full path to directory containing Thermo .raw files 
    assignments_directory : str
        Full path to directory containing CoreMS assignment results 
    eic_tolerance : float
        Tolerance (ppm) for extraction of ion chromatogram (from Thermo .raw file) of given m/z 
    internal_std_mz : float
        m/z of internal standard used for quality control checks; defaults to 678.2915, which is the mass of [cyanocobalamin]2+ (vitamin B12)
    sample_list : str
        Full path to .csv file containing sample list. This file will be imported as Pandas DataFrame 
    [DEPRECATED] csvfile_addend : str
        Textual difference between the name of the Thermo .raw file and the .csv file containing the unprocessed CoreMS assignments. If the .raw file is named SAMPLE.raw, then the assignments file is assumed to be named 'SAMPLE' + csvfile_addend + '.csv'. Defaults to '_assignments'. 
    [DEPRECATED] dispersity_addend : str
        Textual difference between the name of the Thermo .raw file and the .csv file with the processed CoreMS assignments AND dispersity calculation.'
    """

    raw_file_directory: str = ''
    assignments_directory: str = ''
    eic_tolerance: float = 5.0 # ppm 
    internal_std_mz: float = 678.2915 # defaults to mass of [cyanocobalamin]2+
    sample_list: str = '' #  
    time_interval = 2
    std_time_range = [0,20]
    blank_sample_name = ''
    