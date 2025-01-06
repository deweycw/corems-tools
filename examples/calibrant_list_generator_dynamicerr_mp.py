import os

from tempfile import tempdir
import warnings

warnings.filterwarnings('ignore')
from pathlib import Path
import sys
sys.path.append('./')

module_dir = 'c:/Users/15087/Desktop/CoreMSdev/CoreMS'
os.environ['COREMS_DATABASE_URL'] = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
sys.path.append(module_dir)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import multiprocessing as mp
from scipy import stats
import numpy as np

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

def assign_formula(file, time_min, time_max, interval, data_dir, error):
    times = list(range(time_min,time_max,interval))    

    print("\n\nLoading file: "+ file)

    # Read in sample list and load MS data
    MSfiles={}
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)

    MSfiles[file]=parser

    tic=parser.get_tic(ms_type='MS',smooth=False,peak_detection=False)[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})


    MSParameters.molecular_search.min_ppm_error = error - 3 # minimum ppm error for detecting m/z calibrant peaks
    MSParameters.molecular_search.max_ppm_error = error + 3 # maximum ppm error for detecting m/z calibrant peaks

    results = []
    for timestart in times:
        print(timestart)
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()

        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)

        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 20

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 12)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 2)


        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False
        mass_spectrum.molecular_search_settings.max_oc_filter=1.2
        mass_spectrum.molecular_search_settings.max_hc_filter=3
        mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'D': 1,
                                                                        'H': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'Si':4}

        SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)

        assignments=mass_spectrum.to_dataframe()
        assignments['Time']=timestart
        results.append(assignments)

    results=pd.concat(results,ignore_index=True)

    return(results)


def assign_formula(file, times, refmasslist=None): 

	print("Loading file: "+ file)
	
	parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)
	parser.chromatogram_settings.scans = (-1, -1)
	tic=parser.get_tic(ms_type='MS',smooth=False,peak_detection=False)[0]
	tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
	
	results = []
	for t in range(len(times)-1):
		print(times[t])
		t1=times[t]
		t2=times[t+1]
		scans=tic_df[tic_df.time.between(t1,t2)].scan.tolist()
		mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)
		
		if refmasslist:
			MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()
		'''
		mass_spectrum.molecular_search_settings.min_dbe = 0
		mass_spectrum.molecular_search_settings.max_dbe = 20
		mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 40)
		mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 80)
		mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 15)
		mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
		mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 2)
		#mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 1)

		mass_spectrum.molecular_search_settings.isProtonated = True
		mass_spectrum.molecular_search_settings.isRadical = False
		mass_spectrum.molecular_search_settings.isAdduct = False
		mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'H': 1,
                                                                        'D': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'S': 2,
																		'Na': 1
                                                                        }
        '''	
		SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
		mass_spectrum.percentile_assigned(report_error=True)
		
		assignments=mass_spectrum.to_dataframe()
		assignments['Time']=t1
		results.append(assignments)
		
	results=pd.concat(results,ignore_index=True)
	results['file'] = file
	results['Molecular Class']=results['Molecular Formula'].str.replace('\d+', '').str.replace(' ', '')
	results['Molecular Class'][results['Heteroatom Class']=='unassigned']='unassigned'
	results['Molecular Class'][results['Is Isotopologue']==1]='Isotope'
	return(results)


def get_ci(df,fname, ci, it):

    from scipy.optimize import curve_fit

    # pip install uncertainties, if needed
    try:
        import uncertainties.unumpy as unp
        import uncertainties as unc
    except:
        try:
            from pip import main as pipmain
        except:
            from pip._internal import main as pipmain
        pipmain(['install','uncertainties'])
        import uncertainties.unumpy as unp
        import uncertainties as unc

    # import data
    x = df['m/z'].values
    y = df['m/z Error (ppm)'].values
    n = len(y)

    def f(x, a, b, c, d):
        return a*x**3 + b*x**2 + c*x + d
    print(x)
    print(y)
    popt, pcov = curve_fit(f, x, y)

    # retrieve parameter values
    a = popt[0]
    b = popt[1]
    c = popt[2]
    d = popt[3]

    # calculate parameter confidence interval
    a,b,c,d = unc.correlated_values(popt, pcov)

    # calculate regression confidence interval
    px = x
    py = a*px**3 + b*px**2 + c*px + d
    nom = unp.nominal_values(py)
    std = unp.std_devs(py)

    def predband(x, xd, yd, p, func, conf=ci):
        alpha = 1.0 - conf    # significance
        N = xd.size          # data sample size
        var_n = len(p)  # number of parameters
        # Quantile of Student's t distribution for p=(1-alpha/2)
        q = stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
        # Stdev of an individual measurement
        se = np.sqrt(1. / (N - var_n) * \
                    np.sum((yd - func(xd, *p)) ** 2))
        # Auxiliary definitions
        sx = (x - xd.mean()) ** 2
        sxd = np.sum((xd - xd.mean()) ** 2)
        # Predicted values (best-fit model)
        yp = func(x, *p)
        # Prediction band
        dy = q * se * np.sqrt(1.0+ (1.0/N) + (sx/sxd))
        # Upper & lower prediction bands.
        lpb, upb = yp - dy, yp + dy
        return lpb, upb

    lpb, upb = predband(px, x, y, popt, f, conf=ci)

    x_safe = []
    y_safe = []
    x_removed = []
    y_removed = []

    for i in range(len(x)):
        if (y[i] < upb[i]) and (y[i] > lpb[i]):
            x_safe.append(x[i])
            y_safe.append(y[i])
        else:
            x_removed.append(x[i])
            y_removed.append(y[i])

    px = np.arange(min(x),max(x),0.01)
    py = a*px**3 + b*px**2 + c*px + d
    nom = unp.nominal_values(py)
    lpb, upb = predband(px, x, y, popt, f, conf=ci)

    fig, ax = plt.subplots()
    # plot the regression
    ax.scatter(x_safe,y_safe,color = 'C0')
    ax.scatter(x_removed, y_removed, color = 'red')
    ax.plot(px, nom, c='black',marker='o',linewidth=0, markersize = 0.75) #,label='y=a*x^3 + b*x^2 + c*x + d')

    ci_p = ci * 100
    # uncertainty lines (95% confidence)
    ##ax.plot(px, nom - 1.96 * std, c='orange',\
    ##        label= str(int(ci_p)) +'% Confidence Region')
    ##ax.plot(px, nom + 1.96 * std, c='orange')
    # prediction band (95% confidence)
    ax.plot(px, lpb, 'C1o', linewidth=0,markersize=0.5) #,label=str(int(ci_p)) +'% Prediction Band' )
    ax.plot(px, upb, 'C1o', linewidth=0,markersize=0.5)
    ax.set_ylabel('m/z Error (ppm)')
    ax.set_xlabel('m/z')
    #ax.legend(loc='best')

    fig.savefig(fname.split('.')[0] + f'_regression{it}.png', dpi=200)

    filtered_df = df[df['m/z'].isin(x_safe)]
    std_cal = np.std(df['m/z Error (ppm)'])
    return filtered_df, std_cal


def make_cal_list(output, polarity, refmasslist, f, expected_spread, ci=0.68, max_it = 2):

    cal_list = output[output['Confidence Score'] > 0.6]
    cal_list = cal_list[cal_list['Ion Charge'] == 1 * polarity]
    cal_list = cal_list[cal_list['Molecular Class']!='Isotope'].drop_duplicates(subset=['Molecular Formula'])

    it = 1
    cal_list, std_cal = get_ci(cal_list,f, ci, it = it)
    while (std_cal > expected_spread) and (it < max_it):
        it += 1
        cal_list, std_cal = get_ci(cal_list,f, ci, it = it)

    if std_cal > expected_spread:

        warnings.warn('St dev of calibration points exceeds expected spread!')

        return None

    else:
        cal=pd.DataFrame({'# Name':cal_list['Molecular Formula'],
                        'm/z value':cal_list['Calculated m/z'],
                        'charge':cal_list['Ion Charge'],
                        ' ion formula':cal_list['Molecular Formula'],
                        'collision cross section [A^2]':cal_list['Ion Charge']})

        cname = f.replace('.raw','_'+ refmasslist)
        cal.to_csv(cname,sep='\t',index=False)


        fig, ((ax1, ax2)) = plt.subplots(1,2)

        fig.set_size_inches(12, 6)
        sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular Class',data=cal_list,ax=ax1, edgecolor='none')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        ax1.set_title('a', fontweight='bold', loc='left')
        sns.kdeplot(x='m/z Error (ppm)',data=cal_list,hue='Time',ax=ax2,legend=True)
        ax2.set_title('b', fontweight='bold', loc='left')
        fig.tight_layout()
        fname = f.replace('.raw','_calibrants_errorplot.jpg')
        fname = fname.split('/')[-1]

        fig.savefig(fname, dpi=200,format='jpg')

        return cal_list

def make_error_dist_fig(output,f):

    #### Plot and save error distribution figure
    fig, ((ax1, ax2)) = plt.subplots(1,2)
    fig.set_size_inches(12, 6)

    sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular Class',data=output,ax=ax1, edgecolor='none')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax1.set_title('a', fontweight='bold', loc='left')
    sns.kdeplot(x='m/z Error (ppm)',data=output,hue='Time',ax=ax2,legend=True)
    ax2.set_title('b', fontweight='bold', loc='left')
    fig.tight_layout()

    fig.savefig(f.replace('.raw','_errorplot.jpg').split('/')[-1], dpi=200,format='jpg')

# Define CoreMS LCMS functions
def run_calibrants(data_dir,file,refmasslist,error):
    
    # set time intervals
    interval = 4
    time_min = 0
    time_max = 8
    times = list(range(time_min,time_max,interval))
    times.append(time_max)

	# peak picking
    MSParameters.mass_spectrum.min_picking_mz=50	# integer, minimum m/z that will be assigned
    MSParameters.mass_spectrum.max_picking_mz=800 # integer, maximum m/z that will be assigned
    MSParameters.mass_spectrum.noise_threshold_log_nsigma = 10 # integer, should be ~7-10 for most cases (lower will keep more noise)
    MSParameters.ms_peak.peak_min_prominence_percent = 0.05 # relative intensity of lowest vs highest peak in spectrum. Keep around 0.01 to 0.1 . 

	# assigment error and scoring
    MSParameters.molecular_search.min_ppm_error = error - 2 # minimum ppm error for detecting m/z calibrant peaks
    MSParameters.molecular_search.max_ppm_error = error + 2 # maximum ppm error for detecting m/z calibrant peaks   
    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@corems-molformdb-1:5432/coremsapp'
    MSParameters.molecular_search.mz_error_score_weight: float = 0.3
    MSParameters.molecular_search.isotopologue_score_weight: float = 0.7
    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"
	
    # molecular search space
    MSParameters.molecular_search.min_dbe=0
    MSParameters.molecular_search.max_dbe = 20
    MSParameters.molecular_search.usedAtoms['C'] = (1, 50)
    MSParameters.molecular_search.usedAtoms['H'] = (4, 100)
    MSParameters.molecular_search.usedAtoms['O'] = (0, 14)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 2)


    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False
    MSParameters.molecular_search.used_atom_valences = {'C': 4,
                                                        '13C': 4,
                                                        'H': 1,
                                                        'D': 1,
                                                        'O': 2,
                                                        'N': 3,
                                                        }


    output = assign_formula(file = file, times=times)

    #Here, we create a new reference mass list.
    cal_list=output[output['Confidence Score']>.4]
    cal_list=cal_list[cal_list['Ion Charge']==1]
    cal_list=cal_list[cal_list['Is Isotopologue']==0].drop_duplicates(subset=['Molecular Formula'])
        
    #### Plot and save error distribution figure of calibrant list as 'filename_calibrants_errorplot.jpg'
    fig, ((ax1, ax2)) = plt.subplots(1,2)
    fig.set_size_inches(12, 6)
    sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular Class',data=cal_list,ax=ax1, edgecolor='none')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax1.set_title('a', fontweight='bold', loc='left')
    ax1.axhline(y=error, color='black', linestyle='--')  

    sns.kdeplot(x='m/z Error (ppm)',data=cal_list,hue='Time',ax=ax2,legend=True)
    ax2.set_title('b', fontweight='bold', loc='left')
    fig.tight_layout()
    fname = file.replace('.raw','_calibrants_errorplot.jpg')

    fig.savefig(data_dir+fname, dpi=200,format='jpg')

    cal=pd.DataFrame({'# Name':cal_list['Molecular Formula'], 'm/z value':cal_list['Calculated m/z'], 'charge':cal_list['Ion Charge'],' ion formula':cal_list['Molecular Formula'],'collision cross section [A^2]':cal_list['Ion Charge']})

    cname = file.replace('.raw','_'+refmasslist)
    cal.to_csv(data_dir+cname,sep='\t',index=False)

def run_calibrants2(data_dir,file,refmasslist,error):

    # set time intervals
    interval = 2
    time_min = 0
    time_max = 24
    times = list(range(time_min,time_max,interval))
    times.append(time_max)

    polarity = 1

	# peak picking
    MSParameters.mass_spectrum.min_picking_mz=50	# integer, minimum m/z that will be assigned
    MSParameters.mass_spectrum.max_picking_mz=800 # integer, maximum m/z that will be assigned
    MSParameters.mass_spectrum.noise_threshold_log_nsigma = 10 # integer, should be ~7-10 for most cases (lower will keep more noise)
    MSParameters.ms_peak.peak_min_prominence_percent = 0.05 # relative intensity of lowest vs highest peak in spectrum. Keep around 0.01 to 0.1 . 

	# assigment error and scoring
    MSParameters.molecular_search.min_ppm_error = error - 2 # minimum ppm error for detecting m/z calibrant peaks
    MSParameters.molecular_search.max_ppm_error = error + 2 # maximum ppm error for detecting m/z calibrant peaks   
    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@corems-molformdb-1:5432/coremsapp'
    MSParameters.molecular_search.mz_error_score_weight: float = 0.3
    MSParameters.molecular_search.isotopologue_score_weight: float = 0.7
    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"
	
    # molecular search space
    MSParameters.molecular_search.min_dbe=0
    MSParameters.molecular_search.max_dbe = 20
    MSParameters.molecular_search.usedAtoms['C'] = (1, 50)
    MSParameters.molecular_search.usedAtoms['H'] = (4, 100)
    MSParameters.molecular_search.usedAtoms['O'] = (0, 14)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 2)


    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False
    MSParameters.molecular_search.used_atom_valences = {'C': 4,
                                                        '13C': 4,
                                                        'H': 1,
                                                        'D': 1,
                                                        'O': 2,
                                                        'N': 3,
                                                        }


    output = assign_formula(file = file, times=times)

    #Here, we create a new reference mass list.

    cal_list = make_cal_list(output,polarity,refmasslist,file,expected_spread=0.5,ci=0.68,max_it=1)

    if cal_list is not None:

        make_error_dist_fig(output,file)

    else:
        print('No list created for: '+file)


if __name__ == '__main__':

    #Set directories here
    data_dir = 'C:/Users/15087/Desktop/Coremsdev/PT Cocultures/'
    refmasslist = "calibrants.ref"
    samplelist_file = "Coculture_samplelist.csv"
    method=run_calibrants2

    samplelist=pd.read_csv(data_dir+samplelist_file)

    results = []

    os.chdir(data_dir)
    samplelist=pd.read_csv(data_dir+samplelist_file)

    #Run first file alone to ensure that the docker database is generated.
    #for index, row in samplelist.iloc[1:].iterrows():
    #    method(data_dir,samplelist.loc[index,'File'],refmasslist,samplelist.loc[index,'m/z error (ppm)'])
    

    #Run first file alone to ensure that the docker database is generated.
    method(data_dir,samplelist.loc[0,'File'],refmasslist,samplelist.loc[0,'m/z error (ppm)'])

    args=[]
    #for index, row in samplelist.iloc[0:].iterrows():
    for index, row in samplelist.iloc[1:].iterrows():
        file=row.loc['File'] 
        error=row.loc['m/z error (ppm)']
        args.append((data_dir,file,refmasslist,error))
    with mp.Pool(processes=5) as pool:
        results = pool.starmap(method, args)
