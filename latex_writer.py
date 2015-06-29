import os
import time

__author__ = 'mwelland'
__version__ = 1.3
__version_date__ = '11/02/2015'
""" This will be a class to receive a list of objects and a file name
    and compose the output to be written to file
    This will also check if an existing file has the same name and file
    location as the intended output file, and offer to cancel the write
    process or to delete the existing file contents to make way for new
    output
"""

class LatexWriter:

    def __init__(self):
        pass

    @property
    def get_version(self):
        """
        Quick function to grab version details for final printing
        :return:
        """
        return 'Version: {0}, Version Date: {1}'.format(str(__version__), __version_date__)

    def run(self, input_list, filename, write_as_latex):
        self.write_as_latex = write_as_latex
        self.input_list = input_list
        self.filename = filename
        
        if self.write_as_latex:
            self.outfile_name = self.filename+'_'+time.strftime("%d-%m-%Y")+\
                            '_'+time.strftime("%H-%M-%S")+'.tex'
            self.pdfname = self.filename+'_'+time.strftime("%d-%m-%Y")+\
                            '_'+time.strftime("%H-%M-%S")+'.pdf'
        else:
            self.outfile_name = self.filename+'_'+time.strftime("%d-%m-%Y")+\
                            '_'+time.strftime("%H-%M-%S")+ '.txt'

        out = open(self.outfile_name, "w")
        self.fill_output_file(out)
        if self.write_as_latex:
            return self.outfile_name, self.pdfname
        else:
            return self.outfile_name

    def fill_output_file(self, out):

        for line in self.input_list:
            print >> out, line
        print 'File written'