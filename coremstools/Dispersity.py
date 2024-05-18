from numpy import average, nan
from pandas import DataFrame, read_csv

from corems.mass_spectra.input import rawFileReader
from coremstools.Parameters import Settings

from pathlib import Path


def CalculateDispersity(path_to_features):
    """
    Method to calculate dispersity metric. 

    Parameters 
    ----------
    path_to_features : str
        The full path to the feature list or CoreMS assignment output (i.e., the file containing the m/zs to extract) 
    """
    rawfile_dir = Settings.raw_file_directory

    time_interval = Settings.time_interval

    def get_dispersity_rt(row, eics):

        mz = row['m/z']
        time = [row['Time'], row['Time'] + time_interval]
        full_chroma = DataFrame({'EIC':eics[0][mz].eic, 'time':eics[0][mz].time})
        tsub_chroma = full_chroma[full_chroma['time'].between(time[0],time[1])]
        tsub_chroma.sort_values(by='EIC',ascending=False)
        
        tsub_chroma['cumulative'] = tsub_chroma.cumsum()['EIC']/tsub_chroma.sum()['EIC']

        n_points = len(tsub_chroma[tsub_chroma['cumulative']<0.5]+1)
        
        if n_points < 3:
            n_points = 3

        peak_chroma = tsub_chroma.head(n_points)
            
        if peak_chroma['EIC'].sum() > 0:
            d = peak_chroma['time'].std()
            t = average(peak_chroma['time'], weights=peak_chroma['EIC']) 

            return d, t
        else:
            return nan, nan
                
    features_df = read_csv(path_to_features)
    sample_list = list(features_df['File'].unique())
    mz_list = list(features_df['m/z'].unique())

    for sample in sample_list:

        rawfile = rawfile_dir + sample
        parser = rawFileReader.ImportMassSpectraThermoMSFileReader(rawfile)

        eics = parser.get_eics(target_mzs=mz_list, tic_data={}, peak_detection=False, smooth=False)

        features_df['Dispersity: ' + sample], features_df['Retention Time: ' + sample] = zip(*features_df.apply(get_dispersity_rt, eics = eics, axis=1))
    
    features_df.to_csv(path_to_features, index=False)


