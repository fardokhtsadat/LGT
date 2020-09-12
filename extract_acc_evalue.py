import pandas as pd

def modify_string(x):
    with open(x) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    substring = "_alternative"
    for i in range(len(content)):
        if substring in content[i]:
            content[i] = content[i][:content[i].rfind(substring)]
        else:
            content[i] = content[i]
    return content


# following groups queries into related species
def get_species(content):
    species_list = ['BOT', 'CONGO', 'DAR', 'GYROMONAS', 'Hexamita', 'IT1', 'MACHU_PICCU', 'MIS2C', 'PIG', 'SOOS4',
                    'TRIMITUS', 'VLADA7']

    species = {'BOT': [], 'CONGO': [], 'DAR': [], 'GYROMONAS': [], 'Hexamita': [], 'IT1': [], 'MACHU_PICCU': [],
               'MIS2C': [], 'PIG': [], 'SOOS4': [], 'TRIMITUS': [], 'VLADA7': []}

    for i in range(len(content)):
        if content[i].split('_')[0] in species_list:
            element = content[i].split('_')[0]
            species.setdefault(element, [])
            species[element].append(content[i])

        species = {k: v for k, v in species.items() if v}

    return species


#the following function retrives info from pickle dataframe for an orthogroup
def hit_acc(x):
    empty_df = pd.DataFrame()
    for i in x:
        pickle_df = '%s_EUK_df.pkl' % (i)
        df = pd.read_pickle(pickle_df)
        for j in range(len(x[i])):
            empty_df = empty_df.append(df.loc[df['qseqid'] == species[i][j]], ignore_index=True)
    empty_df.to_csv('run_result.csv')

content = modify_string('header_OG0000000')
species = get_species(content)
hit_acc(species)
