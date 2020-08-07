# -*- coding: utf-8 -*-
import argparse
import configparser
from reader import Reader
from latex_writer import LatexWriter
from subprocess import call
from primer_module import Primer
import os
import logging

__author__ = 'mwelland'
__version__ = 2.0
__version_date__ = '06/08/2020'
''' 
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


def check_for_required_folders():
    """
    input: files,
    primers: files,
    output: files,
            tex_files: files,
    requirements,
    """
    top_level_files = os.listdir('.')
    assert all([folder_name in top_level_files for folder_name in
               ['input', 'primers', 'output', 'requirements']])
    assert 'tex_files' in os.listdir('output')


def about():
    return """
    \nMatthew Welland, 8, January 2015
    This is a Python program for creating reference sequences
    To create sequences you will need a computer with LaTex installed\n
    Clicking the "Browse..." button will show you the local directory
    From here choose an LRG file you would like to create a reference for
    To identify the correct LRG, find the gene on http://www.lrg-sequence.org/LRG\n
    This program will produce a .PDF document
    This has been done to prevent any issues with later editing
    The document can be annontated by the use of highlighting\n\n
    If there are any faults during execution or output problems
    please contact matthew.welland@bwnft.nhs.uk, WMRGL, Birmingham\n
    
    ─────────▄──────────────▄
    ────────▌▒█───────────▄▀▒▌
    ────────▌▒▒▀▄───────▄▀▒▒▒▐
    ───────▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐
    ─────▄▄▀▒▒▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐
    ───▄▀▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀██▀▒▌
    ──▐▒▒▒▄▄▄▒▒▒▒▒▒▒▒▒▒▒▒▒▀▄▒▒▌
    ──▌▒▒▐▄█▀▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐
    ─▐▒▒▒▒▒▒▒▒▒▒▒▌██▀▒▒▒▒▒▒▒▒▀▄▌
    ─▌▒▀▄██▄▒▒▒▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌
    ─▌▀▐▄█▄█▌▄▒▀▒▒▒▒▒▒░░░░░░▒▒▒▐
    ▐▒▀▐▀▐▀▒▒▄▄▒▄▒▒▒▒▒░░░░░░▒▒▒▒▌
    ▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒░░░░░░▒▒▒▐
    ─▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌
    ─▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▐
    ──▀▄▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▄▒▒▒▒▌
    ────▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀
    ───▐▀▒▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀
    --──▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀▀
    
    \nSo gene\nSuch reference\nWow
    """


def run_parser(input_file):
    file_type = check_file_type(input_file)
    if file_type == 'gbk':
        from GbkParser import GbkParser
        gbk_reader = GbkParser(input_file, app_settings, args.trim_flanking)
        dictionary = gbk_reader.run()
        parser_details = gbk_reader.get_version
    elif file_type == 'lrg':
        from LrgParser import LrgParser
        lrg_reader = LrgParser(input_file, app_settings, args.trim_flanking)
        dictionary = lrg_reader.run()
        parser_details = lrg_reader.get_version

    else:
        raise Exception("Unrecognised file format: {}".format(file_type))

    primer_list = os.listdir('primers')
    if '{}.csv'.format(dictionary['genename']) in primer_list:
        primer_label = Primer(dictionary, os.getcwd())
        dictionary = primer_label.run()

    parser_details = '{} Parser: {}'.format(
        file_type.upper(),
        parser_details
    )

    os.chdir("output")
    for transcript in dictionary['transcripts']:
        version_details = 'RerenceTypeSetter: {}'.format(get_version())
        list_of_versions = [parser_details, version_details]

        lrg_num = "{}t{}".format(
            input_file.split('.')[0].split('/')[1],
            transcript
        )

        input_reader = Reader(
            dictionary,
            transcript,
            args.write_as_latex,
            list_of_versions,
            args.print_clashes,
            file_type,
            lrg_num,
            args.author,
            app_settings
        )
        input_list, nm = input_reader.run()
        if file_type == 'gbk':
            filename = "{}_{}".format(dictionary['genename'], nm)
        else:
            filename = "{}_{}".format(dictionary['genename'], lrg_num)

        writer = LatexWriter(input_list, filename, args.write_as_latex)
        writer_output = writer.run()
        if args.write_as_latex:
            try:
                call(["pdflatex", "-interaction=batchmode", writer_output[0]])
                clean_up(os.getcwd(), writer_output[1])
                move_files(writer_output[0])
            except:
                logging.error('pdflatex call failed', exc_info=True)


def move_files(latex_file):
    os.rename(latex_file, os.path.join('tex_files', latex_file))


def clean_up(path, pdf_file):
    pdf_split = pdf_file.split('_')
    pwd_files = os.listdir(path)
    pdf_files = [doc for doc in pwd_files if doc.endswith('pdf')]
    for target in pdf_files:        
        if target != 'tex_files':
            target_split = target.split('_')
            if target_split[0:3] == pdf_split[0:3] and target_split[-2:] != pdf_split[-2:]:
                os.remove(os.path.join(path, target))

    targets = [doc for doc in pwd_files if doc.split('.')[-1] not in
               app_settings['DEFAULT']['keep_extensions'].split(',')]
    for target in targets:
        if target != 'tex_files':
            os.remove(os.path.join(path, target))


def check_file_type(file_name):
    """
    identifies LGR (.xml) and GenBank (.gk/gbk) files
    will print an error message and exit the application if other files are specified
    """
    if file_name.endswith('.xml'):
        return 'lrg'
    elif file_name.endswith('.gb') or file_name.endswith('.gbk'):
        return 'gbk'
    else:
        raise Exception('This program only works for GenBank and LRG files')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description='Customise reference sequence settings')
    arg_parser.add_argument('-i', dest='input_file', required=True)
    arg_parser.add_argument('--trim', dest='trim_flanking', action='store_false', default=True)
    arg_parser.add_argument('--clashes', dest='print_clashes', action='store_false', default=True)
    arg_parser.add_argument('--text', dest='write_as_latex', action='store_false', default=True)
    arg_parser.add_argument('--config', dest='settings', default='settings/default_settings.ini')
    arg_parser.add_argument('--author', default="mwelland")
    args = arg_parser.parse_args()

    app_settings = configparser.ConfigParser()
    app_settings.read(args.settings)

    check_for_required_folders()

    keep_extensions = ['pdf', 'tex']
    run_parser(args.input_file)

    print("Process has completed successfully")
