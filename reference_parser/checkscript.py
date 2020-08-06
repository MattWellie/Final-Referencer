import os
from subprocess import call


def clean_up():
    filelist = os.listdir('..')
    for name in filelist:
        os.remove(name)
    os.chdir(os.pardir)
    os.rmdir('testdir')


try:
    import Bio
    from Bio import SeqIO
    print('BioPython is installed')
except ImportError: 
    print('BioPython not Installed')


os.mkdir('testdir')
os.chdir('testdir')
texname = 'QWERTY.tex'

with open(texname, 'w') as outfile:
    outfile.write('\\documentclass[12pt]{article}\n')
    outfile.write('\\begin{document}\n')
    outfile.write('\\end{document}\n')

try:
    call(["pdflatex", "-interaction=batchmode", texname])
    print('PDFLaTex is installed')

except:
    print('PDFLaTex not installed')
    
clean_up()
