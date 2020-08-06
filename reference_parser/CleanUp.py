import os
 
'''
This is a file which allows for the reduction of the 
files present in the reference file output folders, removing
the non-essential files which are created by pdflatex
'''
 
base_contents = os.listdir(os.getcwd())
folders = []
contents = []
keep_extensions = ['pdf', 'tex']


def find(removal_targets):
    for target_folder in removal_targets:
        for group in os.listdir(target_folder):
            removal_files = [
                doc for doc in os.listdir(os.path.join(target_folder, group))
                if doc.split('.')[-1] not in keep_extensions
            ]
            clear(os.path.join(target_folder, group), removal_files)


def clear(path, contents):
    for file_to_remove in contents:
        os.remove(os.path.join(path, file_to_remove))


targets = [folder for folder in base_contents if folder[:4] == 'lrg ']
targets.append('output')
find(targets)
