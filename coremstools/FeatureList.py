from coremstools.Align import Align
from coremstools.Consolidate import Consolidate 
from coremstools.Parameters import Settings
import numpy as np
import pandas as pd
class Features:
    """
    Base class for holding CoreMS features across a dataset. 

    Attributes
    ----------
    feature_list_df : pd.DataFrame or None
        DataFrame holding the aligned and processed features. 
        Each row represents a unique feature, and columns include feature 
        properties, averaged measurements, and intensities across samples.
    sample_list : pd.DataFrame
        DataFrame containing metadata about the samples, primarily a 'File' 
        column listing the raw data file names.

    Methods
    -------
    run_alignment()
        Aligns features across dataset. 
    run_gapfill()
        Runs gapfill. 
    flag_errors()
        Identifies features with potentially significant mass measurement errors based on rolling average and standard error of measured m/z.
    flag_blank_features()
        Calculates a 'blank' flag based on the intensity of a specific blank file compared to the maximum intensity in each feature's spectrum.
    stoichiometric_classification()
        Classifies features based on their stoichiometric ratios of C, H, O, N, P, and S.
    export()
        Writes feature list to .csv file. 
    load_csv()
        Loads feature list from .csv file.
    """
    def __init__(self, sample_list):
        """
        Initializes the Features object.

        Parameters
        ----------
        sample_list : pd.DataFrame
            Pandas DataFrame containing a 'File' column with the name of each 
            .raw file in the dataset (not full path).
        """
                
        self.feature_list_df = None
        self.sample_list = sample_list

    def run_alignment(self):
        """
        Aligns features across all samples in the dataset.
        This method populates `self.feature_list_df` with a DataFrame where each row 
        is a unique feature (e.g., Molecular Formula - Retention Time pair) and columns 
        include averaged parameters (m/z, Confidence Score, etc.) and intensities 
        for that feature in each sample.

        Parameters
        ----------
        The choice between the legacy (dictionary-based) and modern (Pandas-based)
        alignment method is determined by `Settings.alignment_legacy` (default 'False'; use Pandas-based).
        """   
        if Settings.alignment_legacy:
            # Call the legacy (dictionary-based) alignment method
            print('Running legacy alignment...')
            self.feature_list_df = Align.run_legacy(self, self.sample_list)

        else:
            # Call the newer (Pandas-based) alignment method
            self.feature_list_df = Align.run(self, self.sample_list)


    def run_consolidation(self, consolidate_var):
        """
        Consolidates features that are close in m/z and retention time.
        
        This method groups features with indistinguishable m/z values at the same 
        retention times, sums their intensities, and flags 
        features that are redundant. If the feature list 
        (`self.feature_list_df`) doesn't exist, it runs alignment first.

        The method adds/updates the following columns in `self.feature_list_df`:
        - 'consolidated': Indicates whether a feature is part of a consolidated group (1 = yes, 0 = no).
        - 'consolidated flag': A flag indicating whether a feature is the primary representative of a consolidated group (0 = primary, 1 = replaced).
        - 'consolidated id': An identifier for the consolidated group.
        - 'replacement pair': Molecular formula residual elements [{original},{primary replacement}] for features that were consolidated.

        Parameters
        ----------
        consolidate_var : str
            The column name to use for selecting the representative feature 
            from a consolidated group (e.g., 'Confidence Score', 'm/z Error (ppm)').

        The `consolidation_width` and `min_samples` parameters for consolidation 
        are taken from `Settings.consolidation_width` and `Settings.consolidation_min_samples` respectively.
            - consolidation_width : str
                Defines the width for m/z matching. 
                Can be "2sigma", "1sigma", or "fwhm".
            - min_samples : int
                Minimum number of samples ('N Samples' column) for a feature 
                to be considered a primary candidate in consolidation.
        """

        # Ensure feature list exists, run alignment if necessary
        if self.feature_list_df is not None:         
            self.feature_list_df = Consolidate.run(self, consolidate_var, self.feature_list_df)

        else:
            self.run_alignment()
            self.feature_list_df = Consolidate.run(self, consolidate_var, self.feature_list_df)
        

    def flag_errors(self, n_iter=3):
        """
        Flags features with statistically unlikely m/z errors.
        
        This method calculates a rolling average of the 'm/z Error (ppm)' 
        and then computes an 'mz error flag'. The flag is the absolute value 
        of the difference of the rolling average error and the average error 
        of the individual feature divided by 4 times the standard deviation 
        of the m/z error for the feature. 
        This process can be iterated to refine the error assessment.

        The method adds/updates the following columns in `self.feature_list_df`:
        - 'rolling error': The rolling mean of 'm/z Error (ppm)'.
        - 'mz error flag': A numerical flag indicating the likelihood of an error. 
                           Higher values suggest a greater deviation from expected 
                           m/z error. A value greater than 1 indicates a statistically 
                           significant deviation (>99.9% confidence interval).

        Parameters
        ----------
        n_iter : int, optional
            Number of iterations to refine the rolling error and error flag. 
            Defaults to 3.
        """
        print('Running holistic m/z error filter...')

        # Sort by 'Calculated m/z' for the rolling window calculation
        self.feature_list_df = self.feature_list_df.sort_values(by=['Calculated m/z'])

        # Initial calculation of rolling error
        # If 'consolidated flag' exists, calculate rolling error only on primary-consolidated features
        if 'consolidated flag' in self.feature_list_df.columns:
            self.feature_list_df['rolling error'] = self.feature_list_df[self.feature_list_df['consolidated flag'] == 0]['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()
            self.feature_list_df['rolling error'].interpolate(method='linear', inplace=True)

        else:
            self.feature_list_df['rolling error'] = self.feature_list_df['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()

        self.feature_list_df['mz error flag'] = abs(self.feature_list_df['rolling error'] - self.feature_list_df['m/z Error (ppm)']) / (4*self.feature_list_df['m/z Error (ppm)_se'])

        # Iteratively refine the rolling error and m/z error flag
        if n_iter>1:
            for i in range(n_iter-1):
                self.feature_list_df['rolling error'] = np.nan
                if 'consolidated flag' in self.feature_list_df.columns:
                    self.feature_list_df['rolling error'] = self.feature_list_df[(self.feature_list_df['consolidated flag'] == 0) & (self.feature_list_df['mz error flag']<1)]['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()
                    self.feature_list_df['rolling error'].interpolate(method='linear', inplace=True)

                else:
                    self.feature_list_df['rolling error'] = self.feature_list_df[self.feature_list_df['mz error flag']<1]['m/z Error (ppm)'].rolling(window=int(len(self.feature_list_df)/50), center=True, min_periods=0).mean()

                #self.feature_list_df['mz error flag'] = abs(self.feature_list_df['rolling error'] - self.feature_list_df['m/z Error (ppm)']) / (4*self.feature_list_df['m/z Error (ppm)_se']*np.sqrt(self.feature_list_df['N Samples']))
                self.feature_list_df['mz error flag'] = abs(self.feature_list_df['rolling error'] - self.feature_list_df['m/z Error (ppm)']) / (4*self.feature_list_df['m/z Error (ppm)_se'])


    def flag_blank_features(self):
        """
        This function calculates and adds to a feature the following columns:
            'Max Intensity': The maximum signal intensity observed across all samples
            'Max Blank': The maximum signal intensity observed across designated blank files.
            'blank': A ratio (Max Blank / Max Intensity). Higher number means greater presence in blanks.

        It uses `Settings.blank_sample_list` to identify blank sample columns.
        If `self.feature_list_df` is None, it runs alignment first.

        """

        print('Flagging blank features...')

        # Ensure the feature list DataFrame exists; run alignment if not.
        if self.feature_list_df is None:

            self.run_alignment()

        blank_col = []

        # Identify columns corresponding to blank samples
        for blank_sample in Settings.blank_sample_list:

            # Remove file extension if present (e.g., .raw)
            if '.' in blank_sample:
            
                blank_sample = blank_sample.split('.')[0]

            # Find intensity columns ("Intensity:<blank_sample>")
            for col in self.feature_list_df.columns:
                
                if blank_sample in col:

                    blank_col.append(col)

        # Calculate 'Max Intensity' across all intensity columns for each feature
        self.feature_list_df['Max Intensity']=self.feature_list_df.filter(regex='Intensity').max(axis=1)
        # Calculate 'Max Blank' across identified blank intensity columns for each feature
        self.feature_list_df['Max Blank']=self.feature_list_df[blank_col].max(axis=1)
        # Calculate the 'blank' ratio: Max Blank Intensity / Max Overall Intensity
        self.feature_list_df['blank']=self.feature_list_df['Max Blank']/self.feature_list_df['Max Intensity']



    def stoichiometric_classification(self):


        """
        Assigns stoichiometric classifications (e.g., Lipid, Carbohydrate, Peptide)
        to features based on their elemental ratios (O/C, H/C, N/C, P/C) and
        Nominal Oxidation State of Carbon (NOSC).       
        After Rivas et al. (2018) doi: 10.1021/acs.analchem.8b00529

        This method adds/ updates the following columns in `self.feature_list_df`:
        - 'Stoichiometric classification'
        - 'O/C', 'H/C', 'N/C', 'P/C', 'N/P' (elemental ratios)
        - 'NOSC', calculated as: NOSC = 4 - (4*C + H - 3*N - 2*O -2*S + 5*P)/C

        """

        print('Determining stoichiometric classifications...')

        # Initialize the classification column
        self.feature_list_df['Stoichiometric classification']='Unclassified'

        self.feature_list_df['O/C']=self.feature_list_df['O']/self.feature_list_df['C']
        self.feature_list_df['H/C']=self.feature_list_df['H']/self.feature_list_df['C']

        contains_N = True
        contains_P = True
        cols_to_remove = []
        
        # Check for presence of N, P, S columns and add them temporarily if missing for calculation purposes
        if not 'N' in self.feature_list_df.columns:
            self.feature_list_df['N']=0
            cols_to_remove = cols_to_remove + ['N','N/C']
            contains_N = False
        if not 'P' in self.feature_list_df.columns:
            self.feature_list_df['P']=0
            cols_to_remove = cols_to_remove + ['P', 'P/C']
            contains_P = False
        if not 'S' in self.feature_list_df.columns:
            self.feature_list_df['S']=0
            cols_to_remove.append('S')

        if (not contains_N) or (not contains_P):
            cols_to_remove.append('N/P')
        
        # Calculate elemental ratios
        self.feature_list_df['N/C']=self.feature_list_df['N']/self.feature_list_df['C']
        self.feature_list_df['P/C']=self.feature_list_df['P']/self.feature_list_df['C']
        self.feature_list_df['N/P']=self.feature_list_df['N']/self.feature_list_df['P']

        # Calculate Nominal Oxidation State of Carbon (NOSC)
        self.feature_list_df['NOSC'] =  4 -(4*self.feature_list_df['C'] 
                                + self.feature_list_df['H'] 
                                - 3*self.feature_list_df['N'] 
                                - 2*self.feature_list_df['O']
                                - 2*self.feature_list_df['S']
                                +5*self.feature_list_df['P']
                                )/self.feature_list_df['C']

        # Apply rules for 'Lipid' classification
        self.feature_list_df.loc[(self.feature_list_df['O/C']<=0.6) & 
                            (self.feature_list_df['H/C']>=1.32) & 
                            (self.feature_list_df['N/C']<=0.126) &
                            (self.feature_list_df['P/C']<0.35)
                            ,'Stoichiometric classification'] = 'Lipid'

        # Apply rules for 'Phospholipid' classification (overrides Lipid if P > 0)
        self.feature_list_df.loc[(self.feature_list_df['O/C']<=0.6) & 
                            (self.feature_list_df['H/C']>=1.32) & 
                            (self.feature_list_df['N/C']<=0.126) &
                            (self.feature_list_df['P/C']<0.35) &
                            (self.feature_list_df['P']>0)
                            ,'Stoichiometric classification'] = 'Phospholipid'

        # Apply rules for 'A-Sugar' (Amino Sugar) classification
        self.feature_list_df.loc[(self.feature_list_df['O/C']>=0.61) & 
                            (self.feature_list_df['H/C']>=1.45) & 
                            (self.feature_list_df['N/C']>0.07) & 
                            (self.feature_list_df['N/C']<=0.2) & 
                            (self.feature_list_df['P/C']<0.3) & 
                            (self.feature_list_df['O']>=3) &
                            (self.feature_list_df['N']>=1)
                            ,'Stoichiometric classification'] = 'A-Sugar'

        # Apply rules for 'Carbohydrate' classification
        self.feature_list_df.loc[(self.feature_list_df['O/C']>=0.8) & 
                            (self.feature_list_df['H/C']>=1.65) & 
                            (self.feature_list_df['H/C']<2.7) &
                            (self.feature_list_df['O']>=3) &
                            (self.feature_list_df['N']==0)
                            ,'Stoichiometric classification'] = 'Carbohydrate'

        # Apply rules for 'Nucleotide' classification
        self.feature_list_df.loc[(self.feature_list_df['O/C']>=0.5) & 
                            (self.feature_list_df['O/C']<1.7) & 
                            (self.feature_list_df['H/C']>1) & 
                            (self.feature_list_df['H/C']<1.8) &
                            (self.feature_list_df['N/C']>=0.2) & 
                            (self.feature_list_df['N/C']<=0.5) & 
                            (self.feature_list_df['N']>=2) &
                            (self.feature_list_df['P']>=1) &
                            (self.feature_list_df['S']==0) &
                            (self.feature_list_df['Calculated m/z']>305) &
                            (self.feature_list_df['Calculated m/z']<523)
                            ,'Stoichiometric classification'] = 'Nucleotide'

        # Apply rules for 'Phytochemical' classification (Lignin-like compounds)
        self.feature_list_df.loc[(self.feature_list_df['O/C']<=1.15) & 
                            (self.feature_list_df['H/C']<1.32) & 
                            (self.feature_list_df['N/C']<0.126) &
                            (self.feature_list_df['P/C']<=0.2) 
                            ,'Stoichiometric classification'] = 'Phytochemical'

        # Apply rules for 'Organosulfur' classification (overrides others if S > 0)
        self.feature_list_df.loc[(self.feature_list_df['S']>0)
                            ,'Stoichiometric classification'] = 'Organosulfur'

        # Apply rules for 'Peptide' classification (first set of criteria)
        self.feature_list_df.loc[(self.feature_list_df['O/C']>0.12) & 
                            (self.feature_list_df['O/C']<=0.6) & 
                            (self.feature_list_df['H/C']>0.9) & 
                            (self.feature_list_df['H/C']<2.5) & 
                            (self.feature_list_df['N/C']>=0.126) & 
                            (self.feature_list_df['N/C']<=0.7) & 
                            (self.feature_list_df['P/C']<0.17) & 
                            (self.feature_list_df['N']>=1)
                            ,'Stoichiometric classification'] = 'Peptide'
        
        # Apply rules for 'Peptide' classification (second set of criteria, broader O/C and H/C)
        self.feature_list_df.loc[(self.feature_list_df['O/C']>0.6) & 
                            (self.feature_list_df['O/C']<=1) & 
                            (self.feature_list_df['H/C']>1.2) & 
                            (self.feature_list_df['H/C']<2.5) & 
                            (self.feature_list_df['N/C']>=0.2) & 
                            (self.feature_list_df['N/C']<=0.7) & 
                            (self.feature_list_df['P/C']<0.17) & 
                            (self.feature_list_df['N']>=1)
                            ,'Stoichiometric classification'] = 'Peptide'

        # Classify isotopologues specifically
        self.feature_list_df.loc[(self.feature_list_df['Is Isotopologue']>0),'Stoichiometric classification']='Isotoplogue'

        # Remove temporarily added element columns if they were not originally present
        for col in cols_to_remove:
            self.feature_list_df.drop(col, axis = 1, inplace=True)

    def export_csv(self, fname):
        """
        Exports the feature list DataFrame to a CSV file.

        Parameters
        ----------
        fname : str
            The name of the CSV file to be saved. The file will be saved in 
            the directory specified by `Settings.assignments_directory`.
        """
        print('writing to .csv...')
        dir = Settings.assignments_directory
        self.feature_list_df.to_csv(dir + fname, index = False) #, single_file = True, header_first_partition_only = True)
        

    def export_parquet(self):
        """
        Exports the feature list DataFrame to a Parquet file.
        """
        print('writing to .parquet...')
        dir = Settings.assignments_directory
        self.feature_list_df.to_parquet(dir + 'feature_list.parquet', compute =True)

    def load_csv(self, fname):
        """
        Loads a feature list from a .csv file.

        Parameters
        ----------
        fname : str
            The name of the .csv file to load.
        """
        print('Loading feature list from .csv...')
        dir = Settings.assignments_directory
        self.feature_list_df = pd.read_csv(dir + fname, low_memory=False)
    