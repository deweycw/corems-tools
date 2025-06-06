from pandas import unique, concat
from numpy import array, zeros, shape, where, log10, log, sqrt
from tqdm import tqdm
from coremstools.Parameters import Settings
import re

class Consolidate:

    def run(self, consolidate_var, features_df):
        """
        Consolidates features in a DataFrame based on m/z, retention time, and a specified variable.
        Features with indistinguishable m/z values at the same retention times are grouped.
        Their intensities are summed, and a representative feature is chosen based on `consolidate_var`.
        Other features in the group are flagged for removal.

        Args:
            consolidate_var (str): The column name to use for selecting the representative feature 
                                   from a consolidated group (e.g., 'Confidence Score', 'm/z Error (ppm)').
                                   'Confidence Score' is recommended for most cases.
            features_df (pd.DataFrame): DataFrame containing the features to be consolidated.

        Note:
            The `consolidation_width` (e.g., "2sigma", "1sigma", "fwhm") and 
            `min_samples` (minimum 'N Samples' for a primary candidate) parameters 
            are taken from `coremstools.Parameters.Settings`.

        Returns:
            pd.DataFrame: The DataFrame with added/updated columns:
                          - 'consolidated': (int) 1 if part of a consolidated group, 0 otherwise.
                          - 'consolidated flag': (int) 1 if the feature is *not* the chosen representative of its group, 0 if it is.
                          - 'consolidated id': (int) A unique identifier for each consolidated group.
                          - 'replacement pair': (list) For consolidated features, shows the elemental difference from the representative.
        """
        consolidation_width=Settings.consolidation_width 
        min_samples=Settings.consolidation_min_samples
        
        # Initialize new columns in the features DataFrame
        features_df['consolidated'] = 0
        features_df['consolidated flag'] = 0
        features_df['consolidated id'] = 0
        features_df['replacement pair'] = 0

        # Determine the factor for calculating m/z range based on consolidation_width
        # This factor relates the mass resolution to the m/z window (dm)
        if consolidation_width == "2sigma":
            factor = 1 / (sqrt(2 * log(2)))  # Corresponds to ~ +/- 2 standard deviations of a Gaussian peak
        elif consolidation_width == "1sigma":
            factor = 1 / (2 * sqrt(2 * log(2))) # Corresponds to ~ +/- 1 standard deviation
        elif consolidation_width == "fwhm":
            factor = 1 / 2 # Corresponds to Full Width at Half Maximum

        # Identify intensity columns (typically prefixed with 'Intensity')
        intensity_cols = list(features_df.filter(regex='Intensity').columns)

        print('running consolidation...')        
        pbar = tqdm(range(len(features_df.index)))
        
        gf_id = 1 # Initialize a unique ID for each consolidated group

        for ix in pbar:
        
            row = features_df.iloc[ix]
            
            # Process the row only if it hasn't been assigned to a consolidated group yet
            if row['consolidated id'] == 0:
                
                resolution = row['Resolving Power'] 
                mass = row['Calibrated m/z']
                time = row['Time']

                # Calculate the m/z delta (dm) based on mass, resolution, and the chosen factor
                dm = factor * (mass / resolution)
                # Define the m/z range for finding matching features
                mrange = [mass - dm, mass + dm]

                # Find features within the m/z range and at the same retention time
                matches = features_df[(features_df['Calibrated m/z'] > mrange[0]) & (features_df['Calibrated m/z'] < mrange[1]) & (features_df['Time'] == time)]
                
                if(len(matches.index) > 1):
                    
                    # Mark all matched features as consolidated and assign them the same group ID
                    features_df.loc[matches.index,'consolidated'] = 1
                    features_df.loc[matches.index, 'consolidated id'] = gf_id
                    gf_id = gf_id + 1

                    # Sum the intensities of all matched features across all intensity columns
                    matches_sum = matches.filter(regex='Intensity').sum(axis=0)

                    # Update the intensity columns for all matched features with the summed intensities
                    features_df.loc[matches.index, intensity_cols] = matches_sum.to_numpy()
                    
                    # Filter matches to find those that meet the min_samples criteria
                    matches_highn = matches[matches['N Samples'] >= min_samples]

                    # If there are features meeting the min_samples criteria
                    if (len(matches_highn.index) > 0):

                        # Select the 'main' or representative feature based on the consolidate_var
                        if consolidate_var == 'm/z Error (ppm)':
                            # For m/z error, the one with the minimum absolute error is 'main'
                            sub = matches.loc[abs(matches[consolidate_var]) != min(abs(matches[consolidate_var])),'Molecular Formula']
                            main = matches.loc[abs(matches[consolidate_var]) == min(abs(matches_highn[consolidate_var])),'Molecular Formula'].values[0]
                        elif consolidate_var == 'mz error flag':
                            # For mz error flag, the one with the minimum flag value is 'main'
                            main = matches.loc[matches[consolidate_var] != min(matches_highn[consolidate_var]),'Molecular Formula']
                            main = matches.loc[matches[consolidate_var] == min(matches_highn[consolidate_var]),'Molecular Formula'].values[0]
                        else:
                            # For other variables (e.g., Confidence Score), the one with the maximum value is 'main'
                            sub = matches.loc[matches[consolidate_var] != max(matches_highn[consolidate_var]), 'consolidated flag']
                            main = matches.loc[matches[consolidate_var] == max(matches_highn[consolidate_var]),'Molecular Formula'].values[0]
                            
                        # Compare each matched molecule with the 'main' molecule to find differences
                        matches['replacement pair']=matches['Molecular Formula'].apply(lambda row:compare_molecules(row,main))

                        # Flag features that are not the 'main' feature as consolidated (i.e., they will be replaced/represented by 'main')
                        features_df.loc[sub.index, 'consolidated flag'] = 1
                        # Store the comparison result (differences)
                        features_df.loc[matches.index,'replacement pair'] = matches['replacement pair']

                    else:
                        # If no features meet the min_samples criteria, flag all matched features
                        features_df.loc[matches.index, 'consolidated flag'] = 1

        return features_df 
    
    def GapFill_experimental(self, features_ddf):
        
        features_ddf['gapfill'] = False
        features_ddf['gapfill flag'] = False
        intensity_cols = [m for m in features_ddf.columns if '.raw' in m]

        def gapfill(row):

            resolution =  row['Resolving Power'] 
            mass = row['Calibrated m/z']
            time = row['Time']

            mrange = [mass*(1-2/resolution), mass*(1+2/resolution)]

            matches = features_ddf[(features_ddf['Calibrated m/z'] > mrange[0]) & (features_ddf['Calibrated m/z'] < mrange[1]) & (features_ddf['Time'] == time)]

            matches_len = 0
            cs_max = 0
            for part in matches.to_delayed():
                part_len = len(part.compute().index)
                #part_cs = max(part.compute()[gapfill_variable])
                matches_len = matches_len + part_len
                ##if part_cs > cs_max:
                  #  cs_max = part_cs
            if(matches_len > 1):

                row['gapfill'] = True
                
                row[intensity_cols] = max(row[intensity_cols])
                
                if row[gapfill_variable] < max(matches[gapfill_variable]):
                    
                    row['gapfill flag'] = True
            
            return row    
        
        features_ddf_2 = features_ddf.apply(lambda x: gapfill( x,), axis = 1)

        return features_ddf_2
    


    def GapFill_experimental_2(self, results):

        print('performing gap fill')

        intensity_cols = list(results.filter(regex='Intensity').columns)

        cols = results.columns

        results = results[[col for col in cols if 'Intensity' not in col] + [col for col in cols if 'Intensity' in col]]

        results.sort_values(['Time','Calibrated m/z'], inplace=True)

        results['Multiple Peaks w/in Uncertainty'] = None

        results['Gapfill ID'] = None

        results['Gapfill Flag'] = None

        results['Gapfill Molecular Formula'] = None

        holder = []

        for time_step in unique(results['Time']):

            time_step_df = results[results['Time'] == time_step].copy()

            n_rows = len(time_step_df)

            mz_array = array([time_step_df['Calibrated m/z'] for i in range(n_rows)])

            R_array = array([time_step_df['Resolving Power'] for i in range(n_rows)])

            FWHM_array = mz_array / R_array

            mz_error_array = FWHM_array + FWHM_array.T

            mz_diff_array = abs(mz_array - mz_array.T)

            gapfill_inds_array = array(mz_diff_array < mz_error_array)

            n_gapfills_inds = where(gapfill_inds_array)

            n_gapfills_array = zeros(shape(mz_array))

            n_gapfills_array[n_gapfills_inds[0], n_gapfills_inds[1]] = 1

            gapfill_sum = n_gapfills_array.sum(axis=0)

            n_gapfills_vector = array([True if i > 1 else None for i in gapfill_sum ])
            
            time_step_df['Multiple Peaks w/in Uncertainty'] = n_gapfills_vector
            
            time_step_df.sort_values(['Time','Calibrated m/z'], inplace=True)

            offset_diag_rows = array([ i for i in range(0,shape(mz_array)[0]-1)])

            offset_diag_cols = array([ j for j in range(1, shape(mz_array)[0])])

            neighboring_mz_diffs = mz_diff_array[offset_diag_rows, offset_diag_cols]
            
            neighboring_mz_err = mz_error_array[offset_diag_rows, offset_diag_cols]

            residual_diff = neighboring_mz_diffs - neighboring_mz_err

            transition_inds = where(residual_diff > 0 )

            n_true_block = 1
                        
            gapfill_column = zeros((shape(mz_array)[0],1))

            pbar = tqdm(range(len(time_step_df)),ncols=100)

            for ix in pbar:
                
                pbar.set_description_str(desc="Adding gapfill ID for timestep %s" %(time_step) , refresh=True)

                gapfill_bool = time_step_df.iloc[ix].loc['Multiple Peaks w/in Uncertainty']

                if gapfill_bool == True:
                    
                    if log10(n_true_block) < 1:

                        add_string = '.00'
                    
                    elif (log10(n_true_block) >= 1) &(log10(n_true_block) < 2):

                        add_string = '.0'

                    else:

                        add_string = '.'

                    gap_id = float(str(time_step) + add_string + str(n_true_block))
                    
                    gapfill_column[ix,0] = gap_id

                    if ix in transition_inds[0]:

                        n_true_block = n_true_block + 1 

            time_step_df['Gapfill ID'] = gapfill_column
            
            gap_ids_list = [id for id in unique(time_step_df['Gapfill ID']) if id > 0]


            pbar = tqdm(gap_ids_list)
            
            for id in pbar:

                pbar.set_description_str(desc="Adding gapfill flag for timestep %s" %(time_step) , refresh=True)

                id_df = time_step_df[time_step_df['Gapfill ID'] == id].copy()

                id_intensity_sum = id_df.filter(regex='Intensity').sum(axis=0)

                id_df[intensity_cols] = id_intensity_sum

                id_max_confidence_score = max(id_df[gapfill_variable])

                id_df.loc[id_df[gapfill_variable] < id_max_confidence_score,'Gapfill Flag'] = True

                gapfilled_mf = id_df.loc[id_df[gapfill_variable] == id_max_confidence_score,'Molecular Formula']

                id_df.loc[id_df[gapfill_variable] < id_max_confidence_score,'Gapfill Molecular Formula'] = gapfilled_mf.iloc[0]

                time_step_df[time_step_df['Gapfill ID'] == id] = id_df

            holder.append(time_step_df)
        return concat(holder)
        

def compare_molecules(a, b):
  """
  Compares two molecular formulas to find common and unique elements.

  Args:
    a: The first molecular formula string.
    b: The second molecular formula string.

  Returns:
    A tuple containing:
      - core_elements: A dictionary of elements common to both molecules.
      - residual_a: A dictionary of elements unique to molecule a.
      - residual_b: A dictionary of elements unique to molecule b.
  """
  # Parse molecule 'a' to get element counts
  elements_a = {}
  # Regex to find element (e.g., C, 13C, Cl) and its optional count
  for match in re.findall(r"(\d*[A-Z][a-z]*)(\d*)", a):
    element = match[0]
    count = int(match[1]) if match[1] else 1
    elements_a[element] = elements_a.get(element, 0) + count

  elements_b = {}
  for match in re.findall(r"(\d*[A-Z][a-z]*)(\d*)", b):
    element = match[0]
    count = int(match[1]) if match[1] else 1
    elements_b[element] = elements_b.get(element, 0) + count

  core_elements = {}
  residual_a = {}
  residual_b = {}

  # Iterate through elements in molecule 'a'
  for element, count_a in elements_a.items():
    if element in elements_b:
      core_elements[element] = min(count_a, elements_b[element])
      if count_a > elements_b[element]:
        residual_a[element] = count_a - elements_b[element]
    else:
      residual_a[element] = count_a

  # Iterate through elements in molecule 'b'
  for element, count_b in elements_b.items():
    if element in elements_a:
      if count_b > elements_a[element]:
        residual_b[element] = count_b - elements_a[element]
    else:
      residual_b[element] = count_b

  return [residual_a, residual_b]