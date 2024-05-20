from pandas import concat
from re import findall
from tqdm import tqdm

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

def add_mol_class(df):
    '''
    Method to adds molecular class to CoreMS assignments & creates a 'Molecular Class' column. Altered assignment .csvs are rewritten with new columns. 
    '''


    heter = _get_heteroatoms(df)
    molclasses = _get_mol_class(heter)
    df2 = _assign_mol_class(df,molclasses)

    return df2
