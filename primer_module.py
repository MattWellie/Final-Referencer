import re
import os
import csv
import Bio

__author__ = 'Matt'
__version__ = 0.1
__version_date__ = '11/04/2015'


class primer:
    #This is the class for performing primer annotation

    def __init__(self):
        self.dict = {}
        self.primer_dict = {}
        self.basepath = {}
        self.primer_files = []
        self.carry_on = False
        self.csv_reader = ''

    @property
    def get_version(self):
        """
        Quick function to grab version details for final printing
        :return:
        """
        return 'Version: {0}, Version Date: {1}'.format(str(__version__), __version_date__)

    def is_primer_present(self):
        #Quick boolean check to see if there is a file saved with the filename
        genename = self.dict['genename']
        for file in self.primer_files:
            filename = file.split('.')[0]
            if genename.lower() == filename.lower():
                self.carry_on = True
                return filename

    def digest_input(self, filename):
        print filename
        #Extract contents of the CSV
        with open(os.path.join('primers', filename+'.csv')) as csvfile:
            reader = csv.DictReader(csvfile)
            exon = 1
            frag_size = 0
            for row in reader:
                if row['Primer Sequences'] == '':
                    pass
                else:
                    seq = row['Primer Sequences'].upper().strip()
                    if seq == '':
                        pass
                    if row['Exon'] != '':
                        exon = row['Exon']
                    direction = row['Direction']
                    if direction == 'R':
                        seq = self.create_reverse_complement(seq)
                    if row['Fragment Size'] == '':
                        if direction == 'R' and frag_size != '':
                            frag = frag_size
                            frag_size = ''
                        else:
                            frag = ''
                    else:
                        frag = row['Fragment Size']
                        frag_size = frag

                    if row['Primer Batch Numbers'] == '':
                        row['Primer Batch Numbers'] = 'Unavailable'
                    if frag == '':
                        constructed_string = 'Primer %s' % (exon+direction)
                    else:
                        constructed_string = 'Primer %s, Frag size = %s' % (exon+direction, frag)
                    self.search_for_seq(seq, constructed_string)

    def search_for_seq(self, seq, construct):
    
        for transcript in self.dict['transcripts']:
            exonlist = self.dict['transcripts'][transcript]['exons'].keys()
            for exon in exonlist:
                try:
                    match = re.search(r'(?i)%s' % \
                        seq, str(self.dict['transcripts'][transcript]['exons'][exon]['sequence']))
                    if match:
                        #print self.dict['transcripts'][transcript]['exons'][exon]['sequence']
                        self.dict['transcripts'][transcript]['exons'][exon]['sequence'] = \
                                    re.sub(r'(?i)%s' % seq, r'\\pdfcomment[date]{%s}\\hl{%s}' % (construct, match.group()),\
                                    str(self.dict['transcripts'][transcript]['exons'][exon]['sequence']))
                        #print self.dict['transcripts'][transcript]['exons'][exon]['sequence']
                        #this = raw_input()
                except MemoryError:
                    print exon
                    print seq
                    
                    
    def create_reverse_complement(self, string):
        new_list = []
        for x in string:
            if x == 'A':
                new_list.append('T')
            if x == 'C':
                new_list.append('G')
            if x == 'G':
                new_list.append('C')
            if x == 'T':
                new_list.append('A')
        new_string = ''.join(new_list)
        return new_string[::-1]

    def run(self, dictionary, basepath):
        #This is the main method
        self.dict = dictionary
        self.basepath = basepath
        self.primer_files = os.listdir(os.path.join(self.basepath, 'primers'))
        filename = self.is_primer_present()
        if self.carry_on == True:
            #other methods
            self.digest_input(filename)
        return self.dict

