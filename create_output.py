import os

folder_path = '/Users/pgutiea/Downloads/LB'

# Get the names of all files in the folder
file_names = os.listdir(folder_path)
# Get the prefix followed by '_' of all files

prefixes_dict = {}
for file_name in file_names:
    prefix = file_name.split('_')[0]
    if prefix not in prefixes_dict:
        prefixes_dict[prefix] = [file_name]
    else:
        prefixes_dict[prefix].append(file_name)


for prefix, files in prefixes_dict.items():
    new_file = f'{folder_path}/{prefix}_all.txt'
    with open(new_file, 'w') as f:
        for file in files:
            with open(f'{folder_path}/{file}', 'r') as f2:
                for line in f2.readlines():
                    if 'Running' not in line and 'size' not in line:
                        f.write(line.strip() + '\n')
    




# for line in open(f'{folder_path}/{file_names[0]}', 'r').readlines():
#     if 'Running' not in line and 'size' not in line:
#         print(line.strip())