import h5py
import os 
import pandas as pd

def create_h5():

    dir = '/Users/christiandewey/Code/corems-tools/test/testdata/'
    h5f = h5py.File("/Users/christiandewey/Code/corems-tools/test/lcms_dataset.hdf5", "w")

    flist = [f for f in os.listdir(dir) if '.csv' in f]
    print(flist)
    for file in flist:

        group_name = file.split('.')[0]

        #grp = h5f.create_group(file.split('.')[0])

        df = pd.read_csv(dir + file)
        df.drop('file', axis=1, inplace=True)
        df_fields = df.columns

        for field in df_fields:
            
            if df[field].dtype == 'object':
                print(field)
                df.fillna({field: "Unassigned"}, inplace=True)
                #df[field] = df[field].astype(h5py.string_dtype(encoding='utf-8', length=None))
                
        
        for time in df['Time'].unique():

            time_sub = df[df['Time'] == time]

            str_cols ={}

            for field in df_fields:
                if time_sub[field].dtype == 'object':
                    print(field)
                    str_cols[field] = list(time_sub[field])
                    time_sub.drop(field, inplace=True)

            '''for field in df_fields:

                if ('Unnamed' not in field) & ('Time' not in field):
                    
                    if df[field].dtype == 'object':
                        df.fillna({field:"Unassigned"},inplace=True)

                    data_col = list(df[field])

                    if '/' in field:
                        field = field.replace('/','_')'''
            print(time_sub.dtypes)
            dset_name = file + '/' + str(time) #+ '/' + field 
            h5f.create_dataset(dset_name, data = time_sub )

    h5f.close()

        



if __name__ == '__main__':

    create_h5()
