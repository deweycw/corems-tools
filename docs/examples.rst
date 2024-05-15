Examples
===========

.. code-block::

    from coremstools.Parameters import Settings
    from coremstools.Assignments import Assignments 
    import pandas as pd
        
    if __name__ == '__main__':
        
        Settings.raw_file_directory = "/Volumes/IQX-Data/"
        Settings.internal_std_mz = 678.2915

        flist = []
        for f in os.listdir(Settings.raw_file_directory):
            if '.raw' in f:
                flist.append(f)

        df = pd.DataFrame({'File':flist})

        raw_assignments = Assignments(sample_df=df)
        raw_assignments.run_internal_std_qc([10,12])