from pandas import Series, DataFrame
from numpy import inf, nan
import matplotlib.pyplot as plt

from seaborn import histplot

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

        _, axs = plt.subplot_mosaic([['a','b']], figsize=(11,5), constrained_layout=True)
        axs['a'].set(xlabel='Time (min)',ylabel='Intensity',title='Internal Standard EIC = '+str(stdmass) + ' m/z')
        
        print('running QC check ...')
        for file in samplelist['File'].unique():
            #try:
            parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
            parser.chromatogram_settings.eic_tolerance_ppm= Settings.eic_tolerance

            EIC=parser.get_eics(target_mzs=[stdmass],tic_data={},peak_detection=False,smooth=False)
            
            df=DataFrame({'EIC':EIC[0][stdmass].eic,'time':EIC[0][stdmass].time})
            df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]
            area[file]=(sum(df_sub['EIC']))
            rt[file]=(df_sub.time[df_sub.EIC==df_sub.EIC.max()].max())
            axs['a'].plot(df_sub['time'],df_sub['EIC']/1e7,label=file[11:])
            print('  ' + file)
            '''except:
                print('--File not found: ' + file)'''

        #axs['a'].get_legend().remove() #(loc='center left', bbox_to_anchor=(1, 0.5))
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
                samplelist.loc[i,'qc_pass']=1

        print(str(samplelist.qc_pass.sum()) + ' pass of ' + str(len(samplelist)) + ' files (i.e., peak area of standard is <= 2x standard deviation of the mean)')

        peak_stdv=samplelist[samplelist.qc_pass==1].qc_area.std()

        print('std dev of area of standard peak: ' + str(round(peak_stdv/peak_mean*100,1))+'%' )
   
        samplelist.replace([inf, -inf], nan, inplace=True)
        histplot(x='qc_area',data=samplelist,ax=axs['b'])
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
        



