from pandas import DataFrame, Series, read_csv, concat
from seaborn import scatterplot, kdeplot, histplot
import matplotlib.pyplot as plt
from re import findall
from tqdm import tqdm 

from Dispersity import CalculateDispersity
from coremstools.Parameters import Settings

from corems.mass_spectra.input import rawFileReader


class AssignmentCalcs:

    def _get_heteroatoms(df):

        def has_numbers(inputString):
            return any(char.isdigit() for char in inputString)
        
        cols = df.columns
        add = False
        heteros = []
        for c in cols:
            if add:
                if has_numbers(c):
                    continue
                elif c == 'Time':
                    break
                else:
                    heteros.append(c)
            if c == 'O':
                add = True
            elif c == 'Time':
                add = False
        return heteros

    def _get_mol_class(add):

        def get_mol_classes(add, base):
            new = []
            remain = []
            new.append(base)
            for i in range(len(add)):
                new.append(base + add[i])
                new2 = []
                remain = add[i+1:]
                for j in remain:
                    new2 = add[i] + j
                    new.append(base + new2)
            return(new)
        
        base = 'CHO'
        molclasses = []
        for i in range(len(add)):
            e = add[i]
            temp = get_mol_classes(add[i:], base = base)
            base = base+e
            molclasses = molclasses + temp


        output = []
        for x in molclasses:
            if x not in output:
                output.append(x)

        output.append('Unassigned')
        return output


    def _assign_mol_class(all_results, molclasses):

        i = 0 # iterable 
        p = 0 # max len holder
        j = 0 # index of max len mol class

        def get_elements(molclass):            
            elements = [] 
            elements = findall('[A-Z][^A-Z]*', molclass)
            return elements
        

        def get_molclass_subset(included_elements, all_elements, all_results):
            
            tdf = all_results
            for e in all_elements:
                try:
                    tdf[e].fillna(0, inplace = True)
                except:
                    pass

            excluded_elements = [e for e in all_elements if e not in included_elements]
            
            for e in included_elements:
                try:
                    tdf = tdf[tdf[e]>0]
                    for j in excluded_elements:
                        try:
                            tdf = tdf[tdf[j]==0]
                        except:
                            pass
                except:
                    pass

            return tdf
        
        for i in range(0,len(molclasses)):
            mc = molclasses[i]
            if (len(mc) > p) and (mc != 'Unassigned'):
                p = len(mc)
                j = i

        all_elements = get_elements(molclasses[j])

        all_results['ID'] = range(0,len(all_results))

        times = all_results['Time'].unique()

        holder = []
        sizenp = 0
        pbar = tqdm(times)

        for t in pbar:

            time_average = all_results[all_results['Time'] == t]

            sizenp = sizenp + len(time_average)

            for m in molclasses:

                if m != 'Unassigned':
                    elements = get_elements(m)

                    sub = get_molclass_subset(elements, all_elements,time_average[~time_average['Molecular Formula'].isna()]) 
                    sub['Molecular Class'] = m    


                elif m == 'Unassigned':
                    sub = time_average[time_average['Molecular Formula'].isna()] 
                    sub['Molecular Class'] = m

                pbar.set_description_str(desc="Assigning molecular class %s at time %s" % (m,t) , refresh=True)
                holder.append(sub)

        return concat(holder)

    def add_mol_class(self):
        '''
        Method to adds molecular class to CoreMS assignments & creates a 'Molecular Class' column. Altered assignment .csvs are rewritten with new columns. 
        '''

        for f in self.sample_df['File']:

            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            df = read_csv(fpath)
            heter = self._get_heteroatoms(df)
            molclasses = self._get_mol_class(heter)
            df2 = self._assign_mol_class(df,molclasses)
            df2.to_csv(fpath, index = False)


    def ErrorPlot(self, assignments, filename, n_molclass):
        """
        Method to produce plots of assignment error. 

        Parameters 
        ----------
        assignments : DataFrame 
            CoreMS assignments, imported as a CSV file. 
        filename : str
            Name of JPG file to be saved.  
        n_molclass : int
            Specifies number of molecular classes to explicitly represent in error plots. If set to 0, all molecular classes will be explicitly represented. If set to a value greater than 0, the first n_molclass molecular classes, sorted from most abundant to least abundant, will be explicitly represented. 
        """
        
        fig, ((ax1, ax2)) = plt.subplots(1,2)
        fig.set_size_inches(12, 6)
        plot_data = assignments.copy()


        if n_molclass > 0:
            from itertools import islice
            temp_dict = {mc: len(plot_data[plot_data['Molecular Class'] == mc]) for mc in plot_data['Molecular Class'].unique() if mc != 'Unassigned'}
            molclass_num = dict(sorted(temp_dict.items(), key=lambda item: item[1], reverse= True))

            def take(n, iterable):
                if len(iterable) >= n:
                    return list(islice(iterable, n))
                else:
                    return list(islice(iterable,len(iterable)))

            first_n_mcs = take(n_molclass,molclass_num.keys())

            #def top_n_only(mc):
            #    if mc not in first_n_mcs:
            #        mc = 'Other'
            #    return mc

            #plot_data['Molecular Class'] = plot_data['Molecular Class'].transform(top_four_only)
        
            plot_data = plot_data[plot_data['Molecular Class'].isin(first_n_mcs)]

        scatterplot(x='m/z', y='m/z Error (ppm)', hue='Molecular Class', s = 3, data=plot_data, ax=ax1, edgecolor='none')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        ax1.set_title('m/z Error v. Measured m/z', fontweight='bold', loc='center', fontsize='medium')
        kdeplot(x='m/z Error (ppm)', data=assignments, hue='Time', ax=ax2, legend=False)
        #ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        ax2.set_title('m/z Error Distribution, Each Time Window', fontweight='bold', loc='center', fontsize='medium' )
        xpltl = -.05
        ypltl = 1.05
        ax1.text(xpltl, ypltl,'a',
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax1.transAxes, fontweight='bold', fontsize = 12)
        ax2.text(xpltl, ypltl,'b',
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax2.transAxes, fontweight='bold', fontsize = 12)
        
        fig.tight_layout()
        fig.savefig(filename, dpi=200,format='jpg')   


    def RTAssignPlot( self, assignments, filename, n_molclass):
            """
            Method to produce plots of assignments classes across chromatographic separation. 

            Parameters 
            ----------
            assignments : DataFrame 
                CoreMS assignments, imported as a CSV file. 
            filename : str
                Name of JPG file to be saved.    
            """
            plot_data = assignments.copy()


            if n_molclass > 0:
                from itertools import islice
                temp_dict = {mc: len(plot_data[plot_data['Molecular Class'] == mc]) for mc in plot_data['Molecular Class'].unique()}
                molclass_num = dict(sorted(temp_dict.items(), key=lambda item: item[1], reverse= True))

                def take(n, iterable):
                    if len(iterable) >= n:
                        return list(islice(iterable, n))
                    else:
                        return list(islice(iterable,len(iterable)))

                first_n_mcs = take(n_molclass+1,molclass_num.keys())

                def top_n_only(mc):
                    if mc not in first_n_mcs:
                        mc = 'Other'
                    return mc

                plot_data['Molecular Class'] = plot_data['Molecular Class'].transform(top_n_only)
            
            assign_summary=[]
            for time in plot_data['Time'].unique():
                current={}
                current['Time']=time
                mclist = list(plot_data['Molecular Class'].unique())
                if 'Other' in mclist:
                    mclist.remove("Other")
                    mclist.append('Other')
                if 'Unassigned' in mclist:
                    mclist.remove('Unassigned')
                    mclist.append('Unassigned')
                for mol_class in mclist:
                    current[mol_class]=len(plot_data[(plot_data['Molecular Class']==mol_class) & (plot_data['Time']==time)])
                assign_summary.append(current)

            df=DataFrame(assign_summary)
            df=df.sort_values(by='Time')

            df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks')
            plt.legend(bbox_to_anchor=(1.05, 1), title = 'Molecular Class', loc=2, borderaxespad=0.,frameon=False)
            plt.savefig(filename,dpi=200, bbox_inches='tight',format='jpg')


    def StandardQC(samplelist, std_timerange, save_file='internal_std.jpg'):
            
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
        
        print('running QC check ...')
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
                print('  ' + file)
            except:
                print('--File not found: ' + file)

        axs['a'].get_legend().remove() #(loc='center left', bbox_to_anchor=(1, 0.5))
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

        histplot(x='qc_area',data=samplelist,ax=axs['b'])
        axs['b'].set_xlabel('Internal Standard Peak Area')

        xpltl = -.25
        ypltl = 0.98
        axs['a'].text(xpltl, ypltl,'a',
            horizontalalignment='center',
            verticalalignment='center',
            transform = axs['a'].transAxes, fontweight='bold', fontsize = 12)
        axs['b'].text(xpltl, ypltl,'b',
            horizontalalignment='center',
            verticalalignment='center',
            transform = axs['b'].transAxes, fontweight='bold', fontsize = 12)
        
        plt.savefig(data_dir + save_file, dpi=300, bbox_inches = 'tight', format='jpg')

        samplelist.to_csv(Settings.raw_file_directory + Settings.sample_list)



