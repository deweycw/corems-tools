from pandas import DataFrame
from seaborn import scatterplot, kdeplot
import matplotlib.pyplot as plt

class AssignmentError:
    """
    Methods to produce plots of assignment error. 
    """
    def ErrorPlot( self, assignments, filename):
        """
        Method to produce plots of assignment error. 

        Parameters 
        ----------
        assignments : DataFrame 
            CoreMS assignments, imported as a CSV file. 
        filename : str
            Name of JPG file to be saved.    
        """
        
        fig, ((ax1, ax2)) = plt.subplots(1,2)
        fig.set_size_inches(12, 6)
        scatterplot(x='m/z', y='m/z Error (ppm)', hue='Molecular Class', data=assignments, ax=ax1, edgecolor='none')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        ax1.set_title('a', fontweight='bold', loc='left')
        kdeplot(x='m/z Error (ppm)', data=assignments, hue='Time', ax=ax2, legend=False)
        ax2.set_title('b', fontweight='bold', loc='left')
        fig.tight_layout()
        fig.savefig(filename, dpi=200,format='jpg')   


    def RTAssignPlot(self, assignments, filename):
        """
        Method to produce plots of assignments classes across chromatographic separation. 

        Parameters 
        ----------
        assignments : DataFrame 
            CoreMS assignments, imported as a CSV file. 
        filename : str
            Name of JPG file to be saved.    
        """

        assign_summary=[]
        for time in assignments['Time'].unique():
            current={}
            current['Time']=time
            for mol_class in assignments['Molecular Class'].unique():
                current[mol_class]=len(assignments[(assignments['Molecular Class']==mol_class) & (assignments['Time']==time)])
            assign_summary.append(current)

        df=DataFrame(assign_summary)
        df=df.sort_values(by='Time')

        df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        plt.savefig(filename, bbox_inches='tight',format='jpg')
