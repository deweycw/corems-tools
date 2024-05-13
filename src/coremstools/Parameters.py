import os
import dataclasses

@dataclasses.dataclass
class Settings:

    raw_file_directory: str = None
    assignments_directory: str = None
    eic_tolerance: float = 5.0 # ppm 
    internal_std_mz: float = 678.2915 # defaults to mass of [cyanocobalamin]2+
    sample_list: str = None # 
    csvfile_addend = 'assignments'
    dispersity_addend = '_dispersity'
