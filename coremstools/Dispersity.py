from numpy import average, nan
from pandas import DataFrame, read_csv

from corems.mass_spectra.input import rawFileReader
from corems.encapsulation.factory.parameters import LCMSParameters

from coremstools.Parameters import Settings

from tqdm import tqdm

#import memory_profiler

class Dispersity:
    """
    Methods to produce plots of assignment error. 
    """
    @staticmethod
    def CalculateDispersity(sample):
        """
        Method to calculate dispersity metric. 

        Parameters 
        ----------
        sample : str 
            Sample name; corresponds to string in sampe list. 
        time_interval : float
            Time interval overwhich MS scans are averaged in CoreMS assignment.    
        """
        
        LCMSParameters.lc_ms.scans=(-1,-1)
        #addend = Settings.csvfile_addend
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
                t = average(peak_chroma['time']) 

                return d, t
            else:
                return nan, nan

        rawfile = sample.replace('.csv', '.raw')
        parser = rawFileReader.ImportMassSpectraThermoMSFileReader( rawfile)
        
        assignments_file = sample
        assignments = read_csv(assignments_file)

        mzs = list(assignments['m/z'].drop_duplicates())
        eics = parser.get_eics(target_mzs=mzs, tic_data={}, peak_detection=False, smooth=False)
        assignments['Dispersity'], assignments['Retention Time'] = zip(*assignments.apply(get_dispersity_rt, eics = eics, axis=1))
        assignments.to_csv(sample.replace('.csv','_dispersity.csv'), index=False)
