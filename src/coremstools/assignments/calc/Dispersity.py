from numpy import average
from pandas import DataFrame, read_csv

from corems.mass_spectra.input import rawFileReader

from coremstools.Parameters import Settings

class Dispersity:

    def CalculateDispersity(self):

        sample_list = Settings.sample_list
        assignments_dir = Settings.assignments_directory
        rawfile_dir = Settings.raw_file_directory

        def get_dispersity_rt(row, eics):

            mz = row['m/z']
            time = [row['Time'], row['Time'] + self.t_int]
            
            full_chroma = DataFrame({'EIC':eics[0][mz].eic, 'time':eics[0][mz].time})
            
            tsub_chroma = full_chroma[full_chroma['time'].between(time[0],time[1])]
            tsub_chroma.sort_values(by='EIC',ascending=False)
            tsub_chroma['cumulative'] = tsub_chroma.cumsum()['EIC']/tsub_chroma.sum()['EIC']

            n_points = len(tsub_chroma[tsub_chroma['cumulative']<0.5]+1)
            if n_points < 3:
                n_points = 3

            peak_chroma = tsub_chroma.head(n_points)

            d = peak_chroma['time'].std()
            t = average(peak_chroma['time'], weights=peak_chroma['EIC']) 

            return d, t 
        

        for sample in sample_list:

            rawfile = rawfile_dir + sample + '.raw'
            parser = rawFileReader.ImportMassSpectraThermoMSFileReader( rawfile)

            assignments_file = assignments_dir + sample + '.csv'
            assignments = read_csv(assignments_file)

            mzs = assignments['m/z'].drop_duplicates()

            eics = parser.get_eics(target_mzs=mzs, tic_data={}, peak_detection=False, smooth=False)

            assignments['Dispersity'], assignments['Retention Time'] = zip(*assignments.apply(get_dispersity_rt, args=(eics), axis=1))

            assignments.to_csv(assignments_dir + sample + '_dispersity.csv')


