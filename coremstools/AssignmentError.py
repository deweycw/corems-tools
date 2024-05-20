from pandas import DataFrame
from seaborn import scatterplot, kdeplot
import matplotlib.pyplot as plt

class AssignmentError:

    def ErrorPlot( assignments, filename, n_molclass):
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

