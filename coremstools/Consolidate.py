from pandas import unique, concat
from numpy import array, zeros, shape, where, log10, log, sqrt
from tqdm import tqdm

class Consolidate:
    
    def run(self, consolidate_var, features_df, consolidation_width = "2sigma"):
        
        features_df['consolidated'] = 0
        features_df['consolidated flag'] = 0
        features_df['consolidated id'] = 0

        if consolidation_width == "2sigma":
            factor = 1 / (sqrt(2 * log(2)))
        elif consolidation_width == "1sigma":
            factor = 1 / (2 * sqrt(2 * log(2)))
        elif consolidation_width == "fwhm":
            factor = 1 / 2

        intensity_cols = list(features_df.filter(regex='Intensity').columns)

        print('running consolidation...')        
        pbar = tqdm(range(len(features_df.index)))
        
        gf_id = 1

        for ix in pbar:
        
            row = features_df.iloc[ix]

            if row['consolidated id'] == 0:
                
                resolution = row['Resolving Power'] 
                mass = row['Calibrated m/z']
                time = row['Time']

                dm = factor * (1 / resolution)
                mrange = [mass * (1 - dm), mass * (1 + dm)]

                matches = features_df[(features_df['Calibrated m/z'] > mrange[0]) & (features_df['Calibrated m/z'] < mrange[1]) & (features_df['Time'] == time)]
                
                if(len(matches.index) > 1):
                    
                    features_df.loc[matches.index,'consolidated'] = 1
                    features_df.loc[matches.index, 'consolidated id'] = gf_id
                    gf_id = gf_id + 1

                    matches_sum = matches.filter(regex='Intensity').sum(axis=0)

                    features_df.loc[matches.index, intensity_cols] = matches_sum.to_numpy()
                    if consolidate_var == 'm/z Error (ppm)':
                        sub = matches.loc[abs(matches[consolidate_var]) > min(abs(matches[consolidate_var])), 'consolidated flag']
                    elif consolidate_var == 'mz error flag':
                        sub = matches.loc[matches[consolidate_var] > min(matches[consolidate_var]), 'consolidated flag']
                    else:
                        sub = matches.loc[matches[consolidate_var] < max(matches[consolidate_var]), 'consolidated flag']
                    features_df.loc[sub.index, 'consolidated flag'] = 1
        
        return features_df 
    

    def gapfill_legacy(self, featurelist):
        featurelist['gapfill']=False
        featurelist['gapfill flag']=False
        for i, row in featurelist.iterrows():
            resolution=row['Resolving Power']
            mass=row['Calibrated m/z']
            time=row['Time']
            mrange=[mass*(1-2/resolution),mass*(1+2/resolution)]
            matches=featurelist[(featurelist['Calibrated m/z']>mrange[0])&(featurelist['Calibrated m/z']<mrange[1])&(featurelist['Time']==time)]
            if(len(matches)>1):
                featurelist.loc[i,'gapfill']=True
                featurelist.loc[i,featurelist.filter(regex='Intensity').columns]=matches.filter(regex='Intensity').sum(axis=0)
                if featurelist.loc[i,gapfill_variable]<max(matches[gapfill_variable]):
                    featurelist.loc[i,'gapfill flag']=True
        return(featurelist)    

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
        
