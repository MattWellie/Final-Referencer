# -*- coding: utf-8 -*-
import argparse
from Tkinter import *
from tkFileDialog import askopenfilename
from LrgParser import LrgParser
from GbkParser import GbkParser
from reader import Reader
from latex_writer import LatexWriter
from subprocess import call
from primer_module import primer
import os

__author__ = 'mwelland'
__version__ = 1.3
__version_date__ = '11/02/2015'
''' This module of the reference sequence writer creates the user interface.
    This is version 2, for which the individual operational components have
    been abstracted into separate modules.

    Program flow:

    - GUI is generated
        - User chooses an input file (type: LRG (XML) / GenBank
        - User chooses an amount of intronic flanking sequence (number)
        - User clicks 'TRANSLATE'

    - The input file type is checked and the file_type variable is set
        - If the input is LRG, an LRG_Parser instance is created
        - If the input is GenBank, an GbkParser instance is created
        - The appropriate Parser instance is used to read the input file
            contents into a dictionary object which is returned
        - The dictionary has the following structure:

            Dict { pad
                   filename
                   genename
                   refseqname
                   transcripts {  transcript {   protein_seq
                                                 cds_offset
                                                 exons {        exon_number {   genomic_start
                                                                                genomic_stop
                                                                                transcript_start
                                                                                transcript_stop
                                                                                sequence (with pad)

        - Use of this dictionary structure allows for use of absolute references
            to access each required part of the processed input, and allows for
            the extension of the format to include any features required later

    - The returned dictionary is passed through a Reader instance, which scans
        through the created dictionary, and creates a list of Strings which
        represent the typesetting which will be used for the final output.
    - The Reader instance has been chosen to write out in a generic format, to
        allow the dictionary contents to be used as a text output or for LaTex.
        Use of a Boolean write_as_latex variable can be used to decide whether
        the output will include LaTex headers and footers

    - The list output from the Reader instance is written to an output file using
        a writer object. Currently this is a LatexWriter instance, using standard
        printing to file.This could be replaced with a print to .txt for inspection
    - The LatexWriter Class creates an output directory which contains a reference
        to the input file name, intronic padding, the date and time. This is done
        to ensure that the output directory is unique and identifies the exact point
        in time when the output file was created
    - The LatexWriter also creates the full PDF output using a Python facilitated
        command line call. The output '.tex' file is created in the new output
        directory and is processed using pdflatex
'''


def get_version():
    """
    Quick function to grab version details for final printing
    :return:
    """
    return 'Version: {0}, Version Date: {1}'.format(str(__version__), __version_date__)



def open_file():
    current = os.path.join(os.getcwd(), 'input')
    name = askopenfilename(initialdir='%s' % current)
    entry.delete(0, END)
    entry.insert(0, name)


def about():
    print '\nMatthew Welland, 8, January 2015'
    print 'This is a Python program for creating reference sequences'
    print 'To create sequences you will need a computer with LaTex installed\n'
    print 'Clicking the "Browse..." button will show you the local directory'
    print 'From here choose an LRG file you would like to create a reference for'
    print 'To identify the correct LRG, find the gene on http://www.lrg-sequence.org/LRG\n'
    print 'This program will produce a .PDF document'
    print 'This has been done to prevent any issues with later editing'
    print 'The document can be annontated by the use of highlighting\n\n'
    print 'If there are any faults during execution or output problems'
    print 'please contact matthew.welland@bwnft.nhs.uk, WMRGL, Birmingham\n'
    '''
    print '─────────▄──────────────▄'
    print '────────▌▒█───────────▄▀▒▌'
    print '────────▌▒▒▀▄───────▄▀▒▒▒▐'
    print '───────▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐'
    print '─────▄▄▀▒▒▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐'
    print '───▄▀▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀██▀▒▌'
    print '──▐▒▒▒▄▄▄▒▒▒▒▒▒▒▒▒▒▒▒▒▀▄▒▒▌'
    print '──▌▒▒▐▄█▀▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐'
    print '─▐▒▒▒▒▒▒▒▒▒▒▒▌██▀▒▒▒▒▒▒▒▒▀▄▌'
    print '─▌▒▀▄██▄▒▒▒▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌'
    print '─▌▀▐▄█▄█▌▄▒▀▒▒▒▒▒▒░░░░░░▒▒▒▐'
    print '▐▒▀▐▀▐▀▒▒▄▄▒▄▒▒▒▒▒░░░░░░▒▒▒▒▌'
    print '▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒░░░░░░▒▒▒▐'
    print '─▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌'
    print '─▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▐'
    print '──▀▄▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▄▒▒▒▒▌'
    print '────▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀'
    print '───▐▀▒▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀'
    print '--──▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀▀'
    '''
    print '\nSo gene\nSuch reference\nWow'


def run_parser():

    padding = 300
    directory_and_file = entry.get()
    file_name = directory_and_file.split('/')[-2] + '/' + directory_and_file.split('/')[-1]
    file_type = check_file_type(file_name)
    username = entry_name.get()
    dictionary = {}
    nm = ''
    parser_details = ''
    #primer_list = [name.split('.')[0] for name in os.listdir('primers')]
    primer_list = os.listdir('primers')
    primer_applied = False
    if file_type == 'gbk':
        print 'Running parser'
        gbk_reader = GbkParser(file_name, padding, args.trim_flanking)
        dictionary = gbk_reader.run()
        parser_details = gbk_reader.get_version
    elif file_type == 'lrg':
        print 'Running parser'
        lrg_reader = LrgParser(file_name, padding, args.trim_flanking)
        dictionary  = lrg_reader.run()
        parser_details = lrg_reader.get_version

    if dictionary['genename']+'.csv' in primer_list:
        primer_applied = True
        primer_label = primer()
        primer_details = 'Primer Labels: ' + primer_label.get_version
        dictionary = primer_label.run(dictionary, os.getcwd())

    parser_details = '{0} {1} {2}'.format(file_type.upper(), 'Parser:', parser_details)


    os.chdir("output")
    for transcript in dictionary['transcripts']:  
        print 'transcript: %d' % transcript  
        
        input_reader = Reader()
        writer = LatexWriter()
        reader_details = 'Reader: ' + input_reader.get_version
        writer_details = 'Writer: ' + writer.get_version
        xml_gui_details = 'Control: ' + get_version()
        if primer_applied:
            list_of_versions = [parser_details, reader_details, \
                            writer_details, xml_gui_details, primer_details]
        else:
            list_of_versions = [parser_details, reader_details, \
                            writer_details, xml_gui_details]
        lrg_num = file_name.split('.')[0].split('/')[1].replace('_', '\_')+'t'+str(transcript)
        input_list, nm = input_reader.run(dictionary, transcript, args.write_as_latex, list_of_versions, args.print_clashes, file_type, lrg_num, username)
        if file_type == 'gbk':
            filename = dictionary['genename']+'_'+ nm
        else:
            filename = dictionary['genename']+'_'+ file_name.split('.')[0].split('/')[1]+'t'+str(transcript)
        if args.write_as_latex:
            latex_file, pdf_file = writer.run(input_list, filename, args.write_as_latex)
            call(["pdflatex", "-interaction=batchmode", latex_file])
            clean_up(os.getcwd(), pdf_file)
            move_files(latex_file)
        else:
            latex_file = writer.run(input_list, filename, args.write_as_latex)
        
        
        # quick step to allow for non-overlapping writes
        print str(transcript) + ' has been printed'

    print "Process has completed successfully"
    root.quit()
    
def kill_the_spare():
    pass

def move_files(latex):
    os.rename(latex, os.path.join('tex files', latex))
    
def clean_up(path, pdf_file):
    pdf_split = pdf_file.split('_')
    pwd_files = os.listdir(path)
    pdf_files = [doc for doc in pwd_files if \
                doc.split('.')[-1] == 'pdf']
    for target in pdf_files:        
        if target == 'tex files':
            pass  
        else:
            target_split = target.split('_')
            if target_split[0:3] == pdf_split[0:3]\
                and target_split[-2:] != pdf_split[-2:]:
                os.remove(os.path.join(path, target))
    targets = [doc for doc in pwd_files if \
                doc.split('.')[-1] not in keep_extensions]
    for target in targets:
        if target == 'tex files':
            pass   
        else:   
            os.remove(os.path.join(path, target))
        

def check_file_type(file_name):
    """ This function takes the file name which has been selected
        as input. This will identify .xml and .gk/gbk files, and
        will print an error message and exit the application if
        a file is used which does not match either of these types
    """
    if file_name[-4:] == '.xml':
        return 'lrg'
    elif file_name[-3:] == '.gb':
        return 'gbk'
    elif file_name[-4:] == '.gbk':
        return 'gbk'
    else:
        print 'This program only works for GenBank and LRG files'
        exit()

arg_parser = argparse.ArgumentParser(description='Customise reference sequence settings')
arg_parser.add_argument('--trim', dest='trim_flanking', action='store_false', default=True)
arg_parser.add_argument('--clashes', dest='print_clashes', action='store_false', default=True)
arg_parser.add_argument('--text', dest='write_as_latex', action='store_false', default=True)
args=arg_parser.parse_args()

keep_extensions = ['pdf', 'tex']

root = Tk()
menu = Menu(root)
root.config(menu=menu)
helpmenu = Menu(menu)
menu.add_command(label="Help", command=about)

text_in_label = Label(root, text="File name:")
text_in_label.grid(row=0, column=1, sticky='w')
entry = Entry(root)
entry.grid(row=0, column=2, sticky='w')
entry.insert(0, 'input/LRG_292.xml')
button = Button(root, text="Browse...", command=open_file)
button.grid(row=0, column=3)

text2 = Label(root, text="User Name:")
text2.grid(row=3, column=1, sticky='w')
entry_name = Entry(root)
entry_name.grid(row=3, column=2, sticky='w')
entry_name.insert(0, 'Anonymous User')

button = Button(root, text="QUIT", fg="red", command=root.quit)
button.grid(row=4, column=1)
parser = Button(root, text="Translate", fg="blue", command=run_parser)
parser.grid(row=4, column=2)

mainloop()
