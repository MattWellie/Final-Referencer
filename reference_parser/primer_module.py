import re
import os
import csv

from utilities import create_reverse_complement


class Primer:

    def __init__(self, dictionary, basepath):
        self.dict = dictionary
        self.basepath = basepath
        self.primer_dict = {}
        self.primer_files = []
        self.carry_on = False
        self.csv_reader = ''

    def is_primer_present(self):
        # Quick boolean check to see if there is a file saved with the filename
        genename = self.dict['genename']
        for file in self.primer_files:
            filename = file.split('.')[0]
            if genename.lower() == filename.lower():
                self.carry_on = True
                return filename

    def digest_input(self, filename):
        # Extract contents of the CSV
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
                        seq = create_reverse_complement(seq)
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
                    match = re.search(r'(?i){}'.format(seq),
                                      self.dict['transcripts'][transcript]['exons'][exon]['sequence']
                                      )
                    if match:
                        self.dict['transcripts'][transcript]['exons'][exon]['sequence'] = \
                                    re.sub(r'(?i){}'.format(seq), r'\\pdfcomment[date]{%s}\\hl{%s}' % (construct, match.group()),\
                                    str(self.dict['transcripts'][transcript]['exons'][exon]['sequence']))

                except MemoryError:
                    print(exon, seq)

    def run(self):
        self.primer_files = os.listdir(os.path.join(self.basepath, 'primers'))
        filename = self.is_primer_present()
        if self.carry_on:
            self.digest_input(filename)
        return self.dict
