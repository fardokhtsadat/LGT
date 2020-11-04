
import configparser

config = configparser.ConfigParser()
config.read_file(open('config'))

if config['create_pkl_df']['Create_df'] == 'yes':
    wd = config['create_pkl_df']['Working_directory']
    data_path = config['create_pkl_df']['Data_directory'].split(',')
    file_ending = config['create_pkl_df']['File_name_ending'].split(',')

if config['append_existing_pkl_df']['Append'] == 'yes':
    wd = config['append_existing_pkl_df']['Working_directory']
    pkl_df_path = config['append_existing_pkl_df']['Path_to_pkl_df']
    data_path = config['append_existing_pkl_df']['Data_directory'].split(',')
    file_ending = config['append_existing_pkl_df']['File_name_ending'].split(',')

if config['data_selection']['Data_selection'] == 'yes':
    wd = config['data_selection']['Working_directory']
    pkl_df_path = config['data_selection']['Path_to_pkl_df']
    path_input_fasta = config['data_selection']['Path_to_input_file']
    species_names = config['data_selection']['Species_names']
    number_of_top_hits = int(config['data_selection']['Number_of_top_hits'])
    number_of_random_hits = int(config['data_selection']['Number_of_random_hits'])
    db_password = config['data_selection']['DB_password']
    
    



