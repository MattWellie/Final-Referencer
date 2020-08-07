from Bio import SeqIO

__author__ = 'mwelland'
__version__ = 2.0
__version_date__ = '06/08/2020'


class GbkParser:
    """
    Notes:
        Isolated class to deal exclusively with GBK files
        Should return dictionary, not write full output

    Parses the input file to find all the useful values
    This will populate a dictionary to be returned at completion

            Dict { pad
                   genename
                   refseqname
                   transcripts {  transcript {   protein_seq
                                                 cds_offset
                                                 exons {  exon_number {   genomic_start
                                                                          genomic_stop
                                                                          transcript_start
                                                                          transcript_stop
                                                                          sequence (with pad)
    """

    def __init__(self, file_name, app_settings, trim_flanking):

        """
        This class is created by instantiating with a file name and a padding value.
        These are used to locate the appropriate target file, and to select the
        amount of flanking sequence to be appended to exons.
        """
        '''
        :param file_name: the location/identity of the target input file
        :param padding: the required amount of intronic padding
        '''
        padding = int(app_settings['FILE_PARSER']['padding'])
        self.trim_flanking = trim_flanking
        self.exons = []
        self.cds = []
        self.mrna = []
        self.file_name = file_name
        # Read in the specified input file into a variable
        try:
            self.transcriptdict = dict(transcripts={}, input=SeqIO.to_dict(SeqIO.parse(file_name, 'genbank')),
                                       pad=int(padding), pad_offset=int(padding) % 5)
            self.transcriptdict['refseqname'] = self.transcriptdict['input'].keys()[0]
            self.is_matt_awesome = True
        except IOError as fileNotPresent:
            raise Exception("The specified file cannot be located: {}".format(fileNotPresent.filename))

        assert self.transcriptdict['pad'] <= int(app_settings['FILE_PARSER']['max_padding']), \
            "Padding too large, please use a value below {} bases".format(int(app_settings['FILE_PARSER']['padding']))

    @property
    def get_version(self):
        """
        Quick function to grab version details for final printing
        :return:
        """
        return 'Version: {0}, Version Date: {1}'.format(str(__version__), __version_date__)

    def find_cds_delay(self):
        """ Method to find the actual start of the translated sequence
            introduced to sort out non-coding exon problems """
        '''
        :param transcript: currently a relic of the LRG process (29-01-2015), designed to separate the
                           dictionary population process into distinct sections for each transcript
        '''
        for transcript in self.transcriptdict['transcripts'].keys():
            offset_total = 0
            offset = self.transcriptdict['transcripts'][transcript]['cds_offset']
            exon_list = self.transcriptdict['transcripts'][transcript]['list_of_exons']

            for exon in exon_list:
                g_start, g_stop = self.get_genomic_start_end(transcript, exon)
                if offset > g_stop:
                    offset_total = offset_total + (g_stop - g_start)
                elif g_stop > offset >= g_start:
                    self.transcriptdict['transcripts'][transcript]['cds_offset'] = offset_total + (offset - g_start)
                    break

    def get_protein(self):
        """
        This method takes the CDS tagged block from the GenBank features section and parses the
        contents to retrieve the protein sequence. This is added to the appropriate section of
        dictionary used to hold all required details of the file.
        """
        '''
        :param cds: a list containing the cds element(s) of the genbank features
        '''
        for alternative in self.transcriptdict['Alt transcripts']:
            selected_cds = self.cds[alternative-1]
            self.transcriptdict['transcripts'][alternative].update(
                {
                    'protein_seq': '{}*'.format(selected_cds.qualifiers['translation'][0]),
                    'NP_number': selected_cds.qualifiers['protein_id'][0],
                    'cds_offset': selected_cds.location.start
                }
            )
          
    def get_mrna_exons(self):
        """ This uses the list of exon start and stop positions to populate 
            the exon positions in the dictionary""" 

        for alt in self.transcriptdict['Alt transcripts']:
            self.transcriptdict['transcripts'][alt] = {
                'list_of_exons': [],
                'exons': {}
            }
            selected_mrna = self.mrna[alt-1]
            try:
                self.transcriptdict['transcripts'][alt]['NM_number'] = selected_mrna.qualifiers['transcript_id'][0]
            except KeyError:
                self.transcriptdict['transcripts'][alt]['NM_number'] = self.transcriptdict['genename']
                self.transcriptdict['refseqname'] = self.transcriptdict['genename']
                self.transcriptdict['genename'] = self.cds[0].qualifiers['gene'][0]
            exon = 1
            subfeatures = selected_mrna._get_sub_features()
            
            for coords in subfeatures:
                self.transcriptdict['transcripts'][alt]['list_of_exons'].append(exon)
                self.transcriptdict['transcripts'][alt]['exons'][exon] = {
                    'genomic_start': coords.location.start,
                    'genomic_end': coords.location.end
                }
                exon += 1

    def get_exon_contents(self):
        """
        This function is supplied with the list of exon tagged blocks from the features section
        and populates the exons region of the dictionary with the exon number, coordinates, and
        sequence segments which define the exon
        """
        '''
        :param exons: a list of the exon objects from the GenBank features list
        '''
        for alternative in self.transcriptdict['Alt transcripts']:
            sequence = self.transcriptdict['full genomic sequence']
            for exon_number in self.transcriptdict['transcripts'][alternative]['exons'].keys():
                start, end = self.get_genomic_start_end(alternative, exon_number)
                seq = sequence[start:end]
                pad = self.transcriptdict['pad']
                exon_list = self.transcriptdict['transcripts'][alternative]['list_of_exons']
                if pad != 0:
                    if self.trim_flanking:
                        if exon_number < len(exon_list)-1:
                            next_exon = exon_list[exon_number]
                            next_start, _end = self.get_genomic_start_end(alternative, next_exon)
                            if end > next_start-(self.transcriptdict['pad']*2):
                                half_way_point = int(round((next_start - (end+1))/2))
                                if half_way_point % 2 == 1:
                                    half_way_point -= 1
                                pad3 = sequence[end:end+half_way_point]
                            else:
                                assert end + pad <= len(sequence), "Exon index out of bounds"
                                pad3 = sequence[end:end + pad]
                        else:
                            assert end + pad <= len(sequence), "Exon index out of bounds"
                            pad3 = sequence[end:end + pad]

                        if exon_number != exon_list[0]:
                            previous_exon = exon_list[exon_number-2]
                            _start, previous_end = self.get_genomic_start_end(alternative, previous_exon)
                            if start < previous_end+(self.transcriptdict['pad']*2):
                                half_way_point = int(round((start - (previous_end+1))/2))
                                if half_way_point % 2 == 1:
                                    half_way_point -= 1 
                                pad5 = sequence[previous_end+half_way_point+1:start]
                            else:
                                assert start - pad >= 0, "Exon index out of bounds"
                                pad5 = sequence[start - (pad):start]
                        else:
                            assert start - pad >= 0, "Exon index out of bounds"
                            pad5 = sequence[start - (pad):start]
                    else:
                        assert start - pad >= 0, "Exon index out of bounds"
                        assert end + pad <= len(sequence), "Exon index out of bounds"
                        pad3 = sequence[end:end + pad]
                        pad5 = sequence[start - (pad + 1):start - 1]
                        seq = pad5.lower() + seq + pad3.lower()
                    seq = pad5.lower() + seq + pad3.lower()
                self.transcriptdict['transcripts'][alternative]['exons'][exon_number]['sequence'] = seq

    def fill_and_find_features(self):
        dictionary = self.transcriptdict['input'][self.transcriptdict['refseqname']]
        self.transcriptdict['full genomic sequence'] = dictionary.seq
        features = dictionary.features
        for feature in features:
            # Multiple exons are expected, not explicitly used
            if feature.type == 'exon':
                self.exons.append(feature)

        """ 
        This section works on the assumption that each exon in the file will use the appropriate gene name
        and that the only relevant CDS and mRNA sections will also contain the same accession
        """
        try:
            self.transcriptdict['genename'] = self.exons[0].qualifiers['gene'][0]
            for feature in features:
                if feature.type == 'CDS':
                    if feature.qualifiers['gene'][0] == self.transcriptdict['genename']:
                        self.cds.append(feature)
                elif feature.type == 'mRNA':
                    if feature.qualifiers['gene'][0] == self.transcriptdict['genename']:
                        self.mrna.append(feature)
        except KeyError:
            for feature in features:
                if feature.type == 'CDS':
                    self.cds.append(feature)
                elif feature.type == 'mRNA':
                    self.mrna.append(feature)
            note = self.mrna[0].qualifiers['note'][0]
            self.transcriptdict['genename'] = note.split('=')[1]
        assert len(self.cds) == len(self.mrna), "There are a different number of CDS and mRNA"

    def get_genomic_start_end(self, transcript, exon_number):
        """
        this should minimise overall lines in this module
        Args:
            transcript:
            exon_number:

        Returns: start and end values

        """
        start = self.transcriptdict['transcripts'][transcript]['exons'][exon_number]['genomic_start']
        stop = self.transcriptdict['transcripts'][transcript]['exons'][exon_number]['genomic_end']
        return start, stop

    def run(self):
        """
        This is the main method of the GBK Parser. This method is called after class instantiation
        and handles the operation of all the other functions to complete the dictionary which will
        hold all of the sequence and exon details of the gene file being parsed
        """
        '''
        :return transcriptdict: This function fills and returns the dictionary, contents explained in docstring above
        '''
        # initial sequence grabbing and populating dictionaries
        self.fill_and_find_features()
        self.transcriptdict['Alt transcripts'] = range(1, len(self.cds)+1)
        self.get_mrna_exons()
        self.get_protein()
        self.get_exon_contents()
        self.find_cds_delay()
        return self.transcriptdict
