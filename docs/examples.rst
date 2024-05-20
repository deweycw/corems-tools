Examples
===========

CoreMS assignment script
------------------------

.. code-block::

    import os
    import sys
    sys.path.append('./')

    import warnings
    warnings.filterwarnings('ignore')

    from corems.mass_spectra.input import rawFileReader
    from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
    from corems.encapsulation.factory.parameters import MSParameters
    from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

    import pandas as pd

    def assign_formula(file, times, interval): 

        MSParameters.mass_spectrum.min_picking_mz=50
        MSParameters.mass_spectrum.max_picking_mz=800
        MSParameters.ms_peak.peak_min_prominence_percent = 0.02
        MSParameters.molecular_search.min_ppm_error = -1
        MSParameters.molecular_search.max_ppm_error = 1
        MSParameters.mass_spectrum.min_calib_ppm_error = -1
        MSParameters.mass_spectrum.max_calib_ppm_error = 1
        MSParameters.molecular_search_settings.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'     

        parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)
        parser.chromatogram_settings.scans = (-1, -1)

        tic=parser.get_tic(ms_type='MS')[0]
        tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

        ref_file = file.split('.')[0] + '_calibrants_pos.ref'

        MSfiles={}
        MSfiles[file]=parser
        results = []

        for timestart in times:

            print('assiging at ' + str(timestart) + ' min')

            scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
            
            mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)

            MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()

            mass_spectrum.molecular_search_settings.min_dbe = 0
            mass_spectrum.molecular_search_settings.max_dbe = 20

            mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 65)
            mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 88)
            mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 15)
            mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 15)
            mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 1)
            
            mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                            '13C': 4,
                                                                            'H': 1,
                                                                            'D': 1,
                                                                            'O': 2,
                                                                            'N': 3,
                                                                            'S': 2
                                                                            }

            mass_spectrum.molecular_search_settings.isProtonated = True
            mass_spectrum.molecular_search_settings.isRadical = False
            mass_spectrum.molecular_search_settings.isAdduct = False

            SearchMolecularFormulas(mass_spectrum, first_hit=True,ion_charge=1).run_worker_mass_spectrum()
            mass_spectrum.percentile_assigned(report_error=True)
            
            assignments=mass_spectrum.to_dataframe()
            assignments['Time']=timestart

            results.append(assignments)
            
        return(pd.concat(results,ignore_index=True))


    if __name__ == '__main__':

        data_dir = "/Volumes/IQX-Data/"

        results = []

        interval = 2
        time_min = 2
        time_max = 30

        times = list(range(time_min,time_max,interval))

        flist = os.listdir(data_dir)
        f_raw = [f for f in flist if '.raw' in f]
        
        for f in f_raw:
            output = assign_formula(file = data_dir+f, times = times, interval=interval)
            output['file'] = f
            output_name = f.split('.')[0] + '_assignments.csv'
            
            output.to_csv(data_dir+output_name)

CoreMSTools processing script
-----------------------------


.. code-block::
    import warnings
    warnings.filterwarnings('ignore')

    from coremstools.Parameters import Settings
    from coremstools.DataSet import DataSet

    if __name__ == '__main__':
        
        Settings.raw_file_directory = './test/testdata/'
        Settings.assignments_directory = Settings.raw_file_directory
        Settings.internal_std_mz = 678.2915
        Settings.std_time_range = [7,10]
        Settings.time_interval = 2
        Settings.blank_sample_name = '20221103_LBA_Boiteau_Zorbax3p5_qh2o_fullmz'

        dset = DataSet(path_to_sample_list = Settings.raw_file_directory + 'sample_list.csv')

        dset.run_internal_std_qc()

        dset.run_assignment_error_plots()

        dset.run_molclass_retention_plots()

        dset.run_dispersity_calcs()

        dset.run_alignment()

        dset.run_gapfill()

        dset.flag_blank_features()

        dset.export_feature_list()


