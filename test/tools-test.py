import pandas as pd
import os

from coremstools.Parameters import Settings
from coremstools.DataSet import DataSet

if __name__ == '__main__':
    
    Settings.raw_file_directory = "/Volumes/IQX-Data/Oregon State University Data/Dewey/05/test/"
    Settings.internal_std_mz = 678.2915
    flist = []
    for f in os.listdir(Settings.raw_file_directory):
        if '.raw' in f:
            flist.append(f)

    df = pd.DataFrame({'File':flist})
    print(df)

