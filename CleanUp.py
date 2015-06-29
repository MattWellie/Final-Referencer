import os
 
'''
This is a file which allows for the reduction of the 
files present in the reference file output folders, removing
the non-essential files which are created by pdflatex

'''
 
path = os.getcwd()
base_contents = os.listdir(path)
folders = []
contents = []
keep_extensions = ['pdf', 'tex']

def find(targets):
    for folder in targets:
        folders = os.listdir(folder)
        for group in folders:
            contents = [doc for doc in \
                          os.listdir(os.path.join(folder, group))\
                          if doc.split('.')[-1] not in keep_extensions]
            clear(os.path.join(folder, group), contents)

def clear(path, contents):
    print path
    for file in contents:
        os.remove(os.path.join(path, file))
            
targets = [folder for folder in base_contents if folder[:4] == 'lrg ']
targets.append('output')
find(targets)
    