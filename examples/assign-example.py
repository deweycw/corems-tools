import os
import warnings
import pandas as pd

warnings.filterwarnings('ignore')
import sys
sys.path.append('./')

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration


def assign_formula(file, times, interval): 

    MSParameters.molecular_search.min_ppm_error = -0.25
    MSParameters.molecular_search.max_ppm_error = 0.25
    MSParameters.molecular_search.db_chunk_size = 500

    MSParameters.mass_spectrum.min_calib_ppm_error = -1
    MSParameters.mass_spectrum.max_calib_ppm_error = 1
    MSParameters.mass_spectrum.calib_pol_order = 2

    MSParameters.mass_spectrum.min_picking_mz=50
    MSParameters.mass_spectrum.max_picking_mz=800

    MSParameters.ms_peak.peak_min_prominence_percent = 0.02

    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"


    print("Loading file: "+ file)
 
    MSfiles={}
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)
    parser.chromatogram_settings.scans = (-1, -1)
    MSfiles[file]=parser

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    ref_file = file.split('.')[0] + '_calibrants_pos.ref'
    refmasslist = ref_file

    results = []
    for timestart in times:
        print(timestart)
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
        
        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)

        MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()

        mass_spectrum.molecular_search_settings.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp' 

        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 20

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 65)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 88)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 15)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 15)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 1)

        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False
        mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'H': 1,
                                                                        'D': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'S': 2
                                                                        }

    
        SearchMolecularFormulas(mass_spectrum, ion_charge=1).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)
        
        assignments=mass_spectrum.to_dataframe()
        assignments['Time']=timestart
        results.append(assignments)
    
    results=pd.concat(results,ignore_index=True)
    
    return(results)



if __name__ == '__main__':

    data_dir = '/Users/christiandewey/Code/corems-tools/test/testdata/'
    results = []

    interval = 2
    time_min = 10
    time_max = 14
    times = list(range(time_min,time_max,interval))

    flist = os.listdir(data_dir)
    f_raw = [f for f in flist if '.raw' in f]
        
    for f in f_raw:
        output = assign_formula(file = data_dir+f, times = times, interval=interval)
        output['file'] = f
        fname = f.split('.')[0] + '_assignments.csv'
        output.to_csv(data_dir+fname)
