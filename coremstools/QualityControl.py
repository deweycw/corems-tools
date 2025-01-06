from pandas import Series, DataFrame, concat
from numpy import inf, nan
import matplotlib.pyplot as plt
import scipy

from seaborn import histplot
from seaborn import lineplot

from corems.mass_spectra.input import rawFileReader
from corems.encapsulation.factory.parameters import LCMSParameters
from coremstools.Parameters import Settings

class QualityControl:

    def StandardQC(self, samplelist, save_file='internal_std.jpg'):
            
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
        std_timerange = Settings.std_time_range
        LCMSParameters.lc_ms.scans=(-1,-1)
        area={}
        rt={}
        ppmerror={}

        _, axs = plt.subplot_mosaic([['a','b']], figsize=(11,5), constrained_layout=True)
        axs['a'].set(xlabel='Time (min)',ylabel='Intensity',title='Internal Standard EIC = '+str(stdmass) + ' m/z')
        
        print('running QC check ...')
        for file in samplelist['File'].unique():
            parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
            parser.chromatogram_settings.eic_tolerance_ppm= Settings.eic_tolerance

            EIC=parser.get_eics(target_mzs=[stdmass],tic_data={},peak_detection=False,smooth=False)
            
            df=DataFrame({'EIC':EIC[0][stdmass].eic,'time':EIC[0][stdmass].time,'scan':EIC[0][stdmass].scans})
            df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]
            # Baseline estimation: simplest approach is to use the minimum intensity
            baseline = df_sub.EIC.min()            
            area[file]=scipy.integrate.trapz(df_sub['EIC']-baseline,df_sub['time'])
            rt[file]=(df_sub.time[df_sub.EIC==df_sub.EIC.max()].max())
            axs['a'].plot(df_sub['time'],df_sub['EIC']/1e7,label=file[11:])

            #Determine m/z error based on 3 mass spectra from the apex of the internal std peak. 
            scan=df_sub.sort_values('EIC',ascending=False).head(3).scan.to_list()
            mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scan)
            masses=mass_spectrum.to_dataframe()['m/z']
            mdiff=abs(masses-stdmass)
            mass=masses[mdiff==mdiff.min()].to_numpy()
            err=(mass-stdmass)/stdmass*1e6
            ppmerror[file]=err[0]

            print('  ' + file)

        axs['a'].set_ylabel('Intensity (x 1e7)')

        samplelist=samplelist.set_index('File')

        samplelist['QC Area '+str(stdmass)] = Series(area)
        samplelist['QC Retention time '+str(stdmass)] = Series(rt)
        samplelist['m/z error (ppm)'] = Series(ppmerror)

        # Flag outliers with peak area greater than 2x standard deviation of the mean 
        peak_stdv=samplelist['QC Area '+str(stdmass)].std()
        peak_mean=samplelist['QC Area '+str(stdmass)].mean()

        samplelist['QC Pass '+str(stdmass)]=0
        for i in samplelist.index:
            if (abs(samplelist['QC Area '+str(stdmass)][i]-peak_mean)<2*peak_stdv):
                samplelist.loc[i,'QC Pass '+str(stdmass)]=1

        print(str(samplelist['QC Pass '+str(stdmass)].sum()) + ' pass of ' + str(len(samplelist)) + ' files (i.e., peak area of standard is <= 2x standard deviation of the mean)')

        peak_stdv=samplelist[samplelist['QC Pass '+str(stdmass)]==1]['QC Area '+str(stdmass)].std()

        print('std dev of area of standard peak: ' + str(round(peak_stdv/peak_mean*100,1))+'%' )
   
        samplelist.replace([inf, -inf], nan, inplace=True)
        histplot(x='QC Area '+str(stdmass),data=samplelist,ax=axs['b'])
        axs['b'].set_xlabel('Internal Standard Peak Area')
        
        xpltl = -.0
        ypltl = 1.05
        axs['a'].text(xpltl, ypltl,'a',
            horizontalalignment='center',
            verticalalignment='center',
            transform = axs['a'].transAxes, fontweight='bold', fontsize = 12)
        axs['b'].text(xpltl, ypltl,'b',
            horizontalalignment='center',
            verticalalignment='center',
            transform = axs['b'].transAxes, fontweight='bold', fontsize = 12)
        
        plt.savefig(data_dir + save_file, dpi=300, bbox_inches = 'tight', format='jpg')

        samplelist.reset_index(inplace=True)

        return samplelist
    
    def tic_plot(self,samplelist,save_file='TIC_plot.jpg',xlimits=None):
        """
        Plots the total ion chromatogram (TIC) for each sample in the sample list.

        Args:
            samplelist (pandas.DataFrame): A DataFrame containing a 'File' column with file paths.
            filename (str): The filename to save the plot as.

        Returns:
            None
        """
        data_dir = Settings.raw_file_directory

        print('Generating TIC plot ...')

        tics=[]
        for file in samplelist['File'].unique():
            parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
            tic=parser.get_tic(ms_type='MS',peak_detection=False, smooth=False)[0]
            tic_df=DataFrame({'Time': tic.time,'Intensity': tic.tic,'Sample':file.replace('.raw','')})
            tics.append(tic_df)

        tics=concat(tics)
        fig, (ax) = plt.subplots(1)
        lineplot(x='Time',y='Intensity',data=tics,ax=ax, hue='Sample')
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Total Ion Current Intensity')
        if(xlimits):
            ax.set_xlim(xlimits)
        ax.legend()
        plt.tight_layout()
        fig.savefig(data_dir +save_file,dpi=300,format='jpg')