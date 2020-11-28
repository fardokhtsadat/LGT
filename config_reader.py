import configparser, os
  
config = configparser.ConfigParser()
config.read_file(open('config'))

if config['create_pkl_df']['Create_df'] == 'yes':
    wd = config['create_pkl_df']['Working_directory'] #the directory the python scripts are
    print('Your working directory is set to %s' %(wd))
    data_path = config['create_pkl_df']['Data_directory'].split(',')
    print('Blast files are: %s' %(data_path))
    columns = [int(i) for i in config['create_pkl_df']['Columns'].split(',')]
    print('Targeted columns in the blast files are: %s' %(columns))
    output_name = config['create_pkl_df']['Pkl_df_path_and_name']
    print('Output name: %s' %(output_name))
    num_cpu = int(config['create_pkl_df']['Number_of_cpus'])
    print('Number of cpus: %s' %(num_cpu))
    #
    os.chdir(wd)
    from create_pkl_df import *
    print('pkg is imported')
    print('calling the fnc ... ')
    extract_data_main(num_cpu, data_path, output_name, columns)
    print('done!')


if config['append_existing_pkl_df']['Append'] == 'yes':
    wd = config['append_existing_pkl_df']['Working_directory']
    print('Your working directory is set to %s' %(wd))
    pkl_df_path = config['append_existing_pkl_df']['Path_to_pkl_df']
    print('Path to the existing pickled dataframe : %s' %(pkl_df_path))
    data_path = config['append_existing_pkl_df']['Data_directory'].split(',')
    print('Blast files are: %s' %(data_path))
    columns = [int(i) for i in config['append_existing_pkl_df']['Columns'].split(',')]
    print('Targeted columns in the blast files are: %s' %(columns))
    os.chdir(wd)
    from append_df import *
    print('pkg is imported')
    append_pkl_df(pkl_df_path, data_path , columns)

if config['data_selection']['Data_selection'] == 'yes':
    wd = config['data_selection']['Working_directory']
    pkl_df_path = config['data_selection']['Path_to_pkl_df']
    orthogroup = config['data_selection']['Path_to_input_file']
    list_of_names = config['data_selection']['Species_names'].split(',')
    number_of_top_hits = int(config['data_selection']['Number_of_top_hits'])
    number_of_random_hits = int(config['data_selection']['Number_of_random_hits'])
    db_password = config['data_selection']['DB_password']
    output_name = config['data_selection']['Output_path_name']
    os.chdir(wd)
    from LGT import *
    print('pkg is imported')
    wrapper(orthogroup, list_of_names, pkl_df_path, db_password, number_of_top_hits, number_of_random_hits, output_name)


