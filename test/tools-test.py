import pandas as pd
from coremstools.Parameters import Settings
from coremstools.assignments.factory.AssignmentsClass import Assignments 
import os 

if __name__ == '__main__':
    
    Settings.raw_file_directory = "/Volumes/IQX/Oregon State University Data/Dewey/05/test/"
    
    flist = []
    for f in os.listdir(Settings.raw_file_directory):
        if '.raw' in f:
            flist.append(f)

    df = pd.DataFrame({'File':flist})
    print(df)
    raw_assignments = Assignments(sample_df=df)
    raw_assignments.run_internal_std_qc([8,10])