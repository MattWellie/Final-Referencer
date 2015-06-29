try:
    import os
    from subprocess import call
    print 'Python is installed'
except ImportError:
    raise ImportError('This should not be possibe')
    exit()
try:
    import Bio
    from Bio import SeqIO
    print 'BioPython is installed'
except ImportError: 
    raise ImportError('BioPython Not Installed')

def clean_up():
    filelist = os.listdir('.')
    for name in filelist:
        os.remove(name)
    os.chdir(os.pardir)
    os.rmdir('testdir')
    
os.mkdir('testdir')    
os.chdir('testdir')
filename = 'QWERTY'
texname = filename+'.tex'
#print texname
outfile = open(texname, 'w')
print >>outfile, '\\documentclass[12pt]{article}'
print >>outfile, '\\begin{document}'
print >>outfile, '\\end{document}'
outfile.close()

try:
    call(["pdflatex", "-interaction=batchmode", texname])
    print 'PDFLaTex is installed'
except:
    print 'PDFLaTex not installed'
    
clean_up()