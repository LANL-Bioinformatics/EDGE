#######
# Download and install database files
#
#######
import os
import getpass
import pwd

#Prompt user for which files to download, where to download them to and where to uncompress them
tar_files = ['edge_dev_HostIndex.tgz', 'edge_dev_NCBI_genomes.tgz', 'edge_dev_GOTTCHA_db.tgz',
             'edge_dev_amplicons_db.tgz', 'edge_dev_nt_20160426.tgz', 'edge_dev_ShortBRED_Database.tgz',
             'edge_dev_PanGIA_db.tgz', 'edge_dev_diamond_db.tgz'
             ]

def get_download_options():
    total_download_size = 181.325
    acceptable = False
    while not acceptable:
        download_other_host_bwa_index = input('Would you like to download the other Host bwa index (~17.2Gb) for host removal, '
                                    'including pig, sheep, cow, monkey, hamster, and goat? (y/n)')
        if download_other_host_bwa_index == 'y' or download_other_host_bwa_index == 'n':
            acceptable = True
        else:
            print('please enter "y" or "n"!')

    acceptable = False
    while not acceptable:
        download_gottcha2 = input('Would you like to download the GOTTCHA2 databases (23.8Gb, 28.5Gb and 1.7Gb and' 
                                  'contains the custom databases for the GOTTCHA2 taxonomic identification pipeline)? (y/n)')
        if download_gottcha2 == 'y' or download_gottcha2 == 'n':
            acceptable = True
        else:
            print('please enter "y" or "n"!')

    acceptable = False
    while not acceptable:
        download_small_bwa_index = input('Would you like to download the smaller BWA index (this is suggested if your machine has < 32GB'
                                         'of memory)? (y/n)')
        if download_small_bwa_index == 'y' or download_small_bwa_index == 'n':
            acceptable = True
        else:
            print('please enter "y" or "n"!')

    if download_other_host_bwa_index == 'y':
        total_download_size += 51.12
        tar_files.extend(['edge_dev_otherHostIndex.tgz'])

    if download_gottcha2 == 'y':
        total_download_size += 60.2992
        tar_files.extend(['edge_dev_GOTTCHA2_bac_db.tgz', 'edge_dev_GOTTCHA2_euk_db.tgz', 'edge_dev_GOTTCHA2_virus_db.tgz'])

    if download_small_bwa_index == 'y':
        total_download_size -= 18.86
        tar_files.extend(['edge_dev_bwa_mini_index.tgz'])
    else:
        tar_files.extend(['edge_dev_bwa_index.tgz'])

    total_download_size = round(total_download_size,2)

    acceptable = False
    while not acceptable:
        download_directory = input('What directory would you like to download the files to? (total size of download is approximately {0} GB), press ENTER for current directory'.format(total_download_size)) or os.getcwd()
        if os.path.exists(download_directory):
            acceptable = True

    acceptable = False
    while not acceptable:
        database_directory = input('What directory would you like to extract the database files to? Press ENTER for current directory') or os.getcwd()
        if os.path.exists(database_directory):
            acceptable = True

    acceptable = False
    while not acceptable:
        user_name = input('What is your username? (press ENTER to use current user{0})'.format(getpass.getuser())) or getpass.getuser()
        try:
            pwd.getpwnam(user_name)
            acceptable = True
        except:
            print('User does not exist')

    tar_file_string = ','.join(tar_files)
    print('The total download size is {0}\nThe databases being downloaded are: {1}'
          '\nThe .tgz files will be downloaded to {2}\nThe databases will be uncompressed to {3}\n'
          'for user {4}'.format(total_download_size,tar_file_string,download_directory,database_directory, user_name))

    options_ok = input('Are these options ok? (y/n)')
    if options_ok == 'y':
        done = True
    else:
        done = False

    return done

if __name__=='__main__':
    done = False
    while done is False:
        done = get_download_options()