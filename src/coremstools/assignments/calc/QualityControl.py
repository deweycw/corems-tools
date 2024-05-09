from pandas import Series, DataFrame
import matplotlib.pyplot as plt
from seaborn import histplot

from corems.mass_spectra.input import rawFileReader

from coremstools.Parameters import Settings

class QualityControl:

    def StandardQC(self, std_timerange, save_file='internal_std'):
            
        """
        Plots the extracted ion chromatogram (EIC) of the internal standard for each sample,
        calculates the peak area and retention time, flags outliers based on standard deviation,
        and saves the results to a CSV file and plots.

        Args:
            std_timerange (list): A list containing the start and end time (in minutes) of the retention time range for the internal standard peak.
            save_file (str): The filename to save the plot and results as.

        Returns:
            pandas.DataFrame: The sample list with additional columns for QC area, retention time, and QC pass/fail flag.
        """

        data_dir = Settings.raw_file_directory
        stdmass = Settings.internal_std_mz

        area={}
        rt={}

        _, axs = plt.subplot_mosaic([['a','b']], figsize=(11,5), constrained_layout=True)
        axs['a'].set(xlabel='Time (min)',ylabel='Intensity',title='Internal Standard EIC = '+str(stdmass) + ' m/z')

        for file in samplelist['File'].unique():
            try:
                parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
                parser.chromatogram_settings.eic_tolerance_ppm= Settings.eic_tolerance

                EIC=parser.get_eics(target_mzs=[stdmass],tic_data={},peak_detection=False,smooth=False)
                
                df=DataFrame({'EIC':EIC[0][stdmass].eic,'time':EIC[0][stdmass].time})
                df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]
                area[file]=(sum(df_sub['EIC']))
                rt[file]=(df_sub.time[df_sub.EIC==df_sub.EIC.max()].max())
                axs['a'].plot(df_sub['time'],df_sub['EIC']/1e7,label=file[11:])
                print(file)
            except:
                print('No file found: ' + file)

        axs['a'].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axs['a'].set_title('a', fontweight='bold', loc='left')
        axs['a'].set_ylabel('Intensity (x 1e7)')

        samplelist=samplelist.set_index('File')

        samplelist['qc_area'] = Series(area)
        samplelist['QC Retention time'] = Series(rt)

        # Flag outliers with peak area greater than 2x standard deviation of the mean 

        peak_stdv=samplelist.qc_area.std()
        peak_mean=samplelist.qc_area.mean()

        samplelist['qc_pass']=0
        for i in samplelist.index:
            if (abs(samplelist.qc_area[i]-peak_mean)<2*peak_stdv):
                samplelist.qc_pass[i]=1

        print(str(samplelist.qc_pass.sum()) + ' pass of ' + str(len(samplelist)))

        peak_stdv=samplelist[samplelist.qc_pass==1].qc_area.std()

        print(str(round(peak_stdv/peak_mean*100,1))+' % std dev')

        #Create plot of overlaid standard EICs
        histplot(x='qc_area',data=samplelist,ax=axs['b'])
        axs['b'].set_xlabel('Internal Standard Peak Area')
        axs['b'].set_title('b', fontweight='bold', loc='left')

        plt.savefig(save_file,dpi=300,format='jpg')

        return(samplelist)
