
import configparser, os
  
config = configparser.ConfigParser()
config.read_file(open('/home/users/LGT_diplo/config'))

if config['create_pkl_df']['Create_df'] == 'yes':
    wd = config['create_pkl_df']['Working_directory'] #the directory the python scripts are
    data_path = config['create_pkl_df']['Data_directory'].split(',')
    output_name = config['create_pkl_df']['Pkl_df_path_and_name']
    num_cpu = int(config['create_pkl_df']['Number_of_cpus'])
    #
    os.chdir(wd)
    from create_pkl_df import *
    extract_data_main(num_cpu, data_path, output_name)

if config['append_existing_pkl_df']['Append'] == 'yes':
    wd = config['append_existing_pkl_df']['Working_directory']
    pkl_df_path = config['append_existing_pkl_df']['Path_to_pkl_df']
    data_path = config['append_existing_pkl_df']['Data_directory'].split(',')
    file_ending = config['append_existing_pkl_df']['File_name_ending'].split(',')

if config['data_selection']['Data_selection'] == 'yes':
    wd = config['data_selection']['Working_directory']
    pkl_df_path = config['data_selection']['Path_to_pkl_df']
    path_input_fasta = config['data_selection']['Path_to_input_file']
    species_names = config['data_selection']['Species_names'].split(',')
    number_of_top_hits = int(config['data_selection']['Number_of_top_hits'])
    number_of_random_hits = int(config['data_selection']['Number_of_random_hits'])
    db_password = config['data_selection']['DB_password']





