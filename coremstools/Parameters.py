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
    time_interval : int
        Time interval overwhich to average MS scans
    std_time_range  : list
        Min and max time range in which to look for standard; first element of list is min, second element is max
    blank_sample_list : list
        List of blank file names (strings)
    alignment_legacy : bool
        If False (preferred default), uses new alignment method based on concatenated Pandas DataFrame.
        If True, uses legacy alignment method based on dictionary; 
    consolidation_width : str
        Width of the Gaussian used for consolidation, can be "2sigma", "1sigma", or 'fwhm'.
    consolidation_min_samples : int
        Minimum number of samples required to keep a feature in the feature list during consolidation.
    """

    raw_file_directory: str = ''
    assignments_directory: str = ''
    eic_tolerance: float = 5.0 # ppm 
    internal_std_mz: float = 678.2915 # defaults to mass of [cyanocobalamin]2+
    sample_list: str = '' 
    time_interval = 2
    std_time_range = [0,20]
    blank_sample_list = []
    alignment_legacy: bool = False
    consolidation_width: str = "2sigma" # "2sigma", "1sigma", or 'fwhm'
    consolidation_min_samples: int = 1 