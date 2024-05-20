from pandas import DataFrame
import matplotlib.pyplot as plt

class MolClassRetention:

    def RTAssignPlot( self, assignments, jpg_filename, n_molclass):
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
        plt.savefig(jpg_filename,dpi=200, bbox_inches='tight',format='jpg')
