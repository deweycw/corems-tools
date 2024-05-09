from pandas import DataFrame
from seaborn import scatterplot, kdeplot
import matplotlib.pyplot as plt

from coremstools.Parameters import Settings

class AssignmentError:

    def errorplot(self, assignments, filename):
    
        #### Plot and save error distribution figure
        fig, ((ax1, ax2)) = plt.subplots(1,2)
        fig.set_size_inches(12, 6)
        scatterplot(x='m/z', y='m/z Error (ppm)', hue='Molecular Class', data=assignments, ax=ax1, edgecolor='none')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        ax1.set_title('a', fontweight='bold', loc='left')
        kdeplot(x='m/z Error (ppm)', data=assignments, hue='Time', ax=ax2, legend=False)
        ax2.set_title('b', fontweight='bold', loc='left')
        fig.tight_layout()
        fig.savefig(Settings.assignments_directory + filename, dpi=200,format='jpg')   


    def rt_assign_plot(self, LCMS_annotations, filename):

        #### Plot library assignments over time
        assign_summary=[]
        for time in LCMS_annotations['Time'].unique():
            current={}
            current['Time']=time
            for mol_class in LCMS_annotations['Molecular Class'].unique():
                current[mol_class]=len(LCMS_annotations[(LCMS_annotations['Molecular Class']==mol_class) & (LCMS_annotations['Time']==time)])
            assign_summary.append(current)

        df=DataFrame(assign_summary)
        df=df.sort_values(by='Time')

        df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        plt.savefig(Settings.assignments_directory + filename, bbox_inches='tight',format='jpg')