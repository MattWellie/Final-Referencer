import time

__author__ = 'mwelland'
__version__ = 2.0
__version_date__ = '06/08/2020'
""" This will be a class to receive a list of objects and a file name
    and compose the output to be written to file
    This will also check if an existing file has the same name and file
    location as the intended output file, and offer to cancel the write
    process or to delete the existing file contents to make way for new
    output
"""


class LatexWriter:

    def __init__(self, input_list, filename, write_as_latex):
        self.write_as_latex = write_as_latex
        self.input_list = input_list
        self.filename = filename
        self.filename_base = None


    @property
    def get_version(self):
        """
        Quick function to grab version details for final printing
        :return:
        """
        return 'Version: {0}, Version Date: {1}'.format(str(__version__), __version_date__)

    @staticmethod
    def get_date_string():
        return time.strftime("%d-%m-%Y_%H-%M-%S")

    def run(self):
        self.filename_base = "{}_{}".format(
            self.filename,
            self.get_date_string()
        )
        
        if self.write_as_latex:
            self.outfile_name = '{}.tex'.format(self.filename_base)
            self.pdfname = '{}.pdf'.format(self.filename_base)

        else:
            self.outfile_name = '{}.txt'.format(self.filename_base)

        with open(self.outfile_name, "w") as out:
            for line in self.input_list:
                out.write("{}\n".format(line.rstrip()))

        if self.write_as_latex:
            return self.outfile_name, self.pdfname
        else:
            return self.outfile_name
