Workflow
========

Here, we outline a typical workflow for using CoreMS and CoreMSTools to process CoreMS assignments of large, multi-sample LC-MS datasets. The overall goal of the workflow is to generate a feature list. A feature list is a list of assignments that are least likely to be erroneous and, importantly, it includes the abundance of these assignments in each sample in the dataset. In addition to including molecular formulas and the intensities, feature lists also typically include - at a minimum - the measured m/z, assignment error, and retention time of each ion.

The steps of a typical workflow can be divided into three groupings:
    1. Formula assignment with CoreMS. 
    2. Quality control checks
    3. Feature filtering 

Formula assignment with CoreMS
------------------------------

Assignments with CoreMS are performed in Python. It is typically easiest and fastest to write an assignment script and to run the script at the command line. It is also possible to run CoreMS assignments in Jupyter notebook, though performance tends to be slow. 

This section is not intended to provide a detailed description of performing formula assignments with CoreMS. However, we do provide a  some basic example for assigning formulas to a multi-sample LC-MS dataset.

Importing Python modules 
+++++++++++++++++++++++++

.. code-block::

    import os
    import warnings
    import pandas as pd
    import sys

    from corems.mass_spectra.input import rawFileReader
    from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
    from corems.encapsulation.factory.parameters import MSParameters
    from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

    sys.path.append('./')
    warnings.filterwarnings('ignore')

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

