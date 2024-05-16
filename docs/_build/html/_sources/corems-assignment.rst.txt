CoreMS Primer
=============

CoreMSTools requires formula assignments generated with CoreMS. It is typically easiest and fastest to write a standalone Python assignment script, and then to run the script at the command line through a virtual Python environment with all necessary dependencies. It is also possible to run CoreMS assignments in a Jupyter notebook, though performance tends to be slow. 

This section is not intended to provide a detailed description of performing formula assignments with CoreMS, but rather a  basic introduction to assigning formulas on a multi-sample LC-MS dataset, collected on an ultrahigh resolution instrument (FT-ICR-MS, Orbitrap MS).

Importing Python modules 
+++++++++++++++++++++++++

We typically employ a subset of the modules available in CoreMS for our LC-MS analyses. These are listed below. Note that if you are not running a containerized version of CoreMS, you will need to save and run your assignment script from the top-level CoreMS directory. Otherwise Python will not be able to find CoreMS and its external packages.

.. code-block::

    from corems.mass_spectra.input import rawFileReader
    from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
    from corems.encapsulation.factory.parameters import MSParameters
    from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

Defining assignment parameters
++++++++++++++++++++++++++++++

The assignment parameters control and define the set of molecular formulas that could possibly be assigned to individual m/z. As such, these parameters have a profound impact on the assignment results and any interpretation of these results. It is therefore extremely important to select assignment parameters that are reasonable for a given instrument and which reflect the expected (or known) chemistry of the samples. To this end, other chemical information about the samples can be useful for constraining the elements employed in the search. 

Some assignments are set globally (apply to all samples), while others are set for each mass spectrum considered. The list of assignment parameters in CoreMS is extensive, and in many cases their default values do not need to be changed. Here we highlight the parameters that we change most frequently in our analyses. 

Global parameters
~~~~~~~~~~~~~~~~~

Global parameters are set in the MSParameters class. 

.. code-block::

    from corems.encapsulation.factory.parameters import MSParameters

For profile mode data, the m/z range of possible peaks needs to be defined, as does the minimum peak prominence. The minimum peak prominence parameter represents the intensity ratio of the least intense and most intense peaks in a mass spectrum. This roughly corresponds to the expected dynamic range of the instrument. As this value decreases, the number of detected peaks increases. One can easily include far too noise peaks in their assignment routine if this value is set too low. We recommend keeping this value at or above 0.01, though in some cases it may make sense to go lower. 

The default method for calculating noise implements the approach described by `Zhurov et al. (2014) <https://pubs.acs.org/doi/10.1021/ac403278t>`_.

.. code-block::
    
    MSParameters.mass_spectrum.min_picking_mz=50
    MSParameters.mass_spectrum.max_picking_mz=800
    MSParameters.ms_peak.peak_min_prominence_percent = 0.02
 
The range of acceptable assignment error is defined globally by setting minimum and maximum assignment errors. This range should be set to reflect the resolving power and mass accuracy of the instrument. For example, for data collected on the 21T FT-ICR-MS, an acceptable range would be +/- 0.25, whereas for data collected on a high-field Orbitrap, a range of +/-1.5 would be acceptable. Candidate formulas outside of the error range will be rejected. 

.. code-block::

    MSParameters.molecular_search.min_ppm_error = -1
    MSParameters.molecular_search.max_ppm_error = 1

The data should be calibrated before assignments are performed. Calibration of the data corrects for measurement drift across an m/z range. Such drift is typical even for well-calibrated instruments. The default calibration method in CoreMS implements a polynomial correction to the measured m/z using a set of reference m/z. The expected error for the reference masses can be set. Poor instrument calibration can be partially corrected by shifting the ranged of expected mass error for the reference m/z. An acceptable m/z range for the reference m/z is dictated by the instrument.

.. code-block::

    MSParameters.mass_spectrum.min_calib_ppm_error = -1
    MSParameters.mass_spectrum.max_calib_ppm_error = 1

Finally, you will need to define the location of the postgresql database to access or make. If CoreMS is installed locally, use the following. 

.. code-block::

    MSParameters.molecular_search_settings.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp' 
    
If you are running CoreMS through a Docker container, replace 'localhost' with the container name. 

.. code-block::
    
    MSParameters.molecular_search_settings.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@molformdb-1:5432/coremsapp'

Mass spectrum parameters
~~~~~~~~~~~~~~~~~~~~~~~~

Other search parameters are set for each mass spectrum under analysis. To set these parameters, a mass spectrum object needs to be created. Below, we show the creation of a mass spectrum object built by averaging all scans collected between 10 and 12 min. 

.. code-block::

    from corems.mass_spectra.input import rawFileReader

    file = "/Volumes/IQX-Data/my_thermo_data_pos.raw"

    timestart = 10
    interval = 2
    timestop = timestart + interval

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)
    parser.chromatogram_settings.scans = (-1, -1)
    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
    scans=tic_df[tic_df.time.between(timestart,timestop)].scan.tolist()
            
    mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)

This mass spectrum object can then be calibrated against a list of reference m/z. 

.. code-block::

    from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
            
    ref_file = "/Volumes/IQX-Data/mz_refs_pos.raw"

    MzDomainCalibration(mass_spectrum, ref_file, mzsegment=[0,1000]).run()

The search parameters associated with the mass spectrum object generally relate to the formulas that CoreMS should search for in its candidate database (which it creates in Postgres). Often it doesn't make sense to change these parameters between samples that comprise a single dataset, as doing so will produce assignments results that are not strictly comparable across samples. 

Elements included in the search are defined as follows.

.. code-block::

    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 65)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 88)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 15)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 15)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 1)

And the valence to be used in assignments can be defined. CoreMS will account for the valence of all elements included in the search in generating candidate formulas. If you wish to consider multiple valences of a single element, you will need to run the search for each desired valence and rebuild the formula database between assignments. Note that the valence of deuterium (D) must be set explicitly in the current version of CoreMS. 

.. code-block::

    mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                '13C': 4,
                                                                'H': 1,
                                                                'D': 1,
                                                                'O': 2,
                                                                'N': 3,
                                                                'S': 2}

The acceptable range of double-bond equivalents (DBE) can be set as follows. Typically a maximum of 20 DBE is acceptable for small-molecule analysis.  

.. code-block::

    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 20

Finally, the ion type can be defined. Possible ions include (de)protonated species, radical species, and adducts. Typically, we treat all ions as (de)protonated, adducts are defined by inclusion of the adduct in the element list. 

.. code-block::

    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical = False
    mass_spectrum.molecular_search_settings.isAdduct = False


Running the search
++++++++++++++++++

With these parameters set, the search can now be executed on the mass spectrum object. The assignment results can then be exported from the mass spectrum object to a data frame. If you wish to assign multiply charged ions, the `ion_charge` parameter can be changed. Note that this parameter reflects the absolute value of the charge and thus will not differ between positive and negative mode data. 

.. code-block::        

    from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
    
    SearchMolecularFormulas(mass_spectrum, first_hit=True, ion_charge=1).run_worker_mass_spectrum()
    
    assignments=mass_spectrum.to_dataframe()

An example for LC-MS data
+++++++++++++++++++++++++

When assigning LC-MS data, it is necessary to mass spectra corresponding to regular time intervals across the chromatographic separation. An individual mass spectrum object is created by averaging the scans within each interval. To accomplish this, we loop through the time range of the separation at the interval over which we wish to average, creating a mass spectrum object for each time interval, assigning formula to each mass spectrum object, and finally merging the assignments in each time interval into a single dataframe. We do this for each file in our dataset. 

The example below demonstrates how to accomplish this. 

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

    def assign_formula(file, times, interval): # by putting the assignment routine in a function, we can easily loop through and assign the raw files in our dataset

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

        ref_file = file.split('.')[0] + '_calibrants_pos.ref'  # we have a ref_file for each raw file, in the same directory as the raw files

        MSfiles={}
        MSfiles[file]=parser
        results = []

        for timestart in times:  # here we loop through our separation

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
            assignments['Time']=timestart # we add this column to keep track of the time window in which the ions were detected and assigned

            results.append(assignments) # here we add assignments for each time interval to a list
            
        return(pd.concat(results,ignore_index=True)) # here we join the assignments of each time interval into a single dataframe. 


    if __name__ == '__main__':

        data_dir = "/Volumes/IQX-Data/"  # this directory contains all the .raw files in our dataset 

        results = []

        interval = 2
        time_min = 2
        time_max = 30

        times = list(range(time_min,time_max,interval))

        flist = os.listdir(data_dir)
        f_raw = [f for f in flist if '.raw' in f]   # this creates a list of .raw files which we loop through
        
        for f in f_raw:

            output = assign_formula(file = data_dir+f, times = times, interval=interval)
            output['file'] = f
            output_name = f.split('.')[0] + '_assignments.csv'
            
            output.to_csv(data_dir+output_name)  # we save assignment output for each raw fle as a csv