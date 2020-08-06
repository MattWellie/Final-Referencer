import re


''' 
This is the Reader class which uses the completed dictionary
from the Parser classes as input and create a list object
containing each of the lines to be used as output in the
appropriate order to be printed.

Partitioning the output process will allow the program to be
able to output to a simple .txt format or to a LaTex document
with appropriate headers and formatting
This will also allow the output to be formally tested or a
variety of assertions to be performed before attempting to
generate output
'''


class Reader:
    """
    This class is used to create list object which will contain the lines
    to be printed to a final output file

    This separates the creation of the dictionary, the writing of the
    output and the actual creation of the output file
    """

    def __init__(self, dictionary, transcript, write_as_latex, list_of_versions,
                 print_clashes, file_type, filename, username):
        self.username = username
        self.list_of_versions = list_of_versions
        self.transcriptdict = dictionary
        self.filename = filename
        self.write_as_LaTex = write_as_latex
        self.transcript = transcript
        self.print_clashes = print_clashes
        self.file_type = file_type
        self.nm = ''
        self.output_list = []
        self.amino_printing = False
        self.amino_spacing = False
        self.exon_spacing = False
        self.exon_printed = False
        self.dont_print = False
        self.check_AA = True
        self.line_break_print = False
        self.pattern = re.compile(r'\\p.*?l{')
        self.found_first_slash = False
        
        # This is a codon-AA dictionary construction created by Peter Collingridge
        # http://www.petercollingridge.co.uk/book/export/html/474
        bases = ['T', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        self.codon_table = dict(zip(codons, amino_acids))

    def decide_number_string_character(self, char, wait_value, cds_count, amino_acid_counter, post_protein_printer,
                                       intron_offset, intron_in_padding, protein_length, intron_out):
        """
        :param char: the next base character to be written
        :param wait_value: variable indicating whether a number has already been written
        :param cds_count: the CDS position of the current char/base
        :param amino_acid_counter: position of the amino acid in the protein sequence
        :param post_protein_printer: to indicate the position of the exon after stop codon
        :param intron_offset:
        :param intron_in_padding:
        :param protein_length: length of protein sequence
        :param intron_out:
        :return: all input values complete with appropriate additions and subtractions, and
                 the new character(s) to be added to the number string

        This function is responsible for determining the next value(s) to be added to the
        output string which will appear above the base sequence.
        """

        output = ''
        if char.isupper():
            # Upper case
            if amino_acid_counter < protein_length:
                if cds_count % 10 == 1 and wait_value == 0:
                    output = '|' + str(cds_count)
                    wait_value = len(str(cds_count))
                    cds_count += 1
                elif wait_value != 0:
                    cds_count += 1
                    wait_value -= 1
                elif wait_value == 0:
                    output = ' '
                    cds_count += 1
            elif amino_acid_counter >= protein_length:

                if post_protein_printer % 10 == 0 and wait_value == 0:
                    output = '|*' + str(post_protein_printer + 1)
                    wait_value = len(str(post_protein_printer + 1)) + 1
                    post_protein_printer += 1
                elif wait_value != 0:
                    wait_value -= 1
                    post_protein_printer += 1
                elif wait_value == 0:
                    post_protein_printer += 1
                    output = ' '
        else:
            if not self.exon_printed:

                # Intron before exon
                self.dont_print = True
                if intron_offset != 0:
                    intron_offset -= 1
                    intron_in_padding -= 1
                    output = ' '
                elif intron_offset == 0 and intron_in_padding % 5 == 0:
                    output = '.'
                    intron_in_padding -= 1
                elif intron_offset == 0 and intron_in_padding % 5 != 0:
                    output = ' '
                    intron_in_padding -= 1
            elif self.exon_printed:
                # Intron after exon
                if wait_value != 0:
                    wait_value -= 1
                    intron_out += 1
                elif wait_value == 0:
                    if intron_out % 5 == 4:
                        output = '.'
                        intron_out += 1
                    elif intron_out % 5 != 4:
                        output = ' '
                        intron_out += 1

        return (output, wait_value, cds_count, amino_acid_counter,
                post_protein_printer, intron_offset, intron_in_padding,intron_out)

    def decide_amino_string_character(self, char, codon_count, amino_acid_counter, codon_numbered, protein):
        """
        :param char: the next base character to be written
        :param codon_count: the position of the current base within the reading frame
        :param amino_acid_counter: the position of the current AA in the protein sequence
        :param codon_numbered: Boolean; has the current codon been numbered
        :param protein: the full protein sequence
        :return: all input values complete with appropriate additions and subtractions, and
                 the new character(s) to be added to the amino number string

        This function determines the next character to be added to the string which shows
        the numbering for the amino acid labelling line
        """

        output = ''
        if char.isupper() and amino_acid_counter < len(protein):  # Bug, condition added
            if self.amino_printing:
                self.amino_spacing = True
                if codon_count == 3:
                    output = protein[amino_acid_counter:amino_acid_counter + 1][0]
                    amino_acid_counter += 1
                    codon_numbered = False
                    codon_count = 1
                else:
                    codon_count += 1
                    output = ' '
            elif not self.amino_printing:
                output = ' '
        elif char.islower():
            output = ' '

        return output, codon_count, amino_acid_counter, codon_numbered

    def print_latex_header(self, refseqid):
        """
        :param refseqid: reference sequence identifier for current input

        This function can be called if the file to be written out is designed to be
        executed as a LaTex script. This will insert the appropriate preamble to
        allow an article class document to be produced which uses a verbatim output
        operation
        """
        try:
            self.nm = self.transcriptdict['transcripts'][self.transcript]['NM_number']
            rep_nm = self.transcriptdict['transcripts'][self.transcript]['NM_number'].replace('_', '\_')  # Required for LaTex
            np = self.transcriptdict['transcripts'][self.transcript]['NP_number'].replace('_', '\_')  # Required for LaTex
        except KeyError:
            print('Additional details not present')
        
        self.output_list.extend(
            [
                '\\documentclass{article}',
                '\\usepackage{color, soul}',
                '\\usepackage{alltt}',
                '\\usepackage{pdfcomment}'
            ]
        )
        self.print_pdfinfo()
        self.output_list.extend(
            [
                '\\begin{document}',
                '\\begin{center}',
                '\\begin{large}'
            ]
        )

        self.output_list.extend(
            [
                'Gene: {} - Sequence: {}\\\\'.format(
                    self.transcriptdict['genename'],
                    refseqid
                ),
                'Transcript: {} - Protein: {}\n'.format(rep_nm, np)
            ]
        )
        if self.file_type == 'lrg':
            self.output_list.append('LRG: {} - Date : \\today'.format(self.filename))
        else:
            self.output_list.append('Date : \\today')
        self.print_pdfinfo()  
        self.output_list.extend(
            [
                '\\end{large}',
                '\\end{center}',
                '$1^{st}$ line: Base numbering. Full stops for intronic +/- 5, 10, 15...\\\\',
                '$2^{nd}$ line: Base sequence. lower case Introns, upper case Exons\\\\',
                '$3^{rd}$ line: Amino acid sequence. Printed on FIRST base of codon\\\\',
                '$4^{th}$ line: Amino acid numbering. Numbered on $1^{st}$ and increments of 10\\\\',
                '\\begin{alltt}'
            ]
        )

    def print_pdfinfo(self):
        self.output_list.extend(
            [
                '\\hypersetup{pdfauthor={},'.format(self.username),
                'pdftitle={Reference sequence for gene: %s}}' % self.nm
            ]
        )

    def print_latex(self):

        """
        This large class imports the entire dictionary and scans through the input dictionary
        to build the output. This output will be kept as a list of strings (with the appropriate
        order maintained) to allow a later decision on whether to print to LaTex, txt or other.

        This function makes a range of calls out to other functions to determine the numbers and
        formatting to use as it scans through the exon(s) base by base
        """

        latex_dict = self.transcriptdict['transcripts'][self.transcript]
        ''' Creates a LaTex file which can be converted to a final document
            Lengths of numbers calculated using len(#)'''
        protein = latex_dict['protein_seq']
        refseqid = self.transcriptdict['refseqname'].replace('_', '\_')  # Required for LaTex
        assert isinstance(latex_dict, dict)
        # A variable to keep a count of the
        # transcript length across all exons
        cds_count = 1 - latex_dict['cds_offset']

        # Account for the number 0 being skipped
        cds_count -= 1
        # The CDS begins at one, the preceeding base is -1. There is no 0
        # The writer must skip 0, so the extra length compensates to keep
        # values in the correct places

        # The initial line(s) of the LaTex file, required to execute
        if self.write_as_LaTex:
            self.print_latex_header(refseqid)

        lines_on_page = 10
        wait_value = 0
        codon_count = 3  # Print AA at start of codon
        amino_acid_counter = 0  # Begin at AA index 0 (first)
        amino_wait = 0  # No number string printed, no wait yet
        codon_numbered = False  # First AA has not been numbered already
        post_protein_printer = 0  # The number for 3' intron '+###' counting
        exon_list = latex_dict['list_of_exons']
        for position in range(len(exon_list)):
            exon_number = int(latex_dict['list_of_exons'][position])
            intron_offset = self.transcriptdict['pad_offset']
            intron_in_padding = self.transcriptdict['pad']
            intron_out = 0  # Or 0?
            self.exon_printed = False
            number_string = []
            dna_string = []
            amino_string = []
            amino_number_string = []
            self.amino_spacing = False
            self.exon_spacing = False
            exon_dict = latex_dict['exons'][exon_number]
            ex_start = exon_dict['genomic_start']
            ex_end = exon_dict['genomic_end']
            if self.file_type == 'gbk':
                ex_start += 1
            self.output_list.append('Exon {} | Start: {} | End: {} | Length: {}'.format(
                exon_number,
                ex_start,
                ex_end,
                ex_end - ex_start
            ))

            if self.print_clashes:
                """ This section allows for a note to be written where the 'intronic' flanking sequence
                    of an exon contains part of the next exon. This serves to clarify whether any overlap
                    may take place. This will not impact the printed output

                    A companion segment in the parser classes is responsible for altering the flanking region
                    if the regions are to avoid overlaps.
                """
                clash_after = False
                clash_before = False
                if exon_number < len(latex_dict['list_of_exons'])-1:
                    try:
                        if ex_end > latex_dict['exons'][latex_dict['list_of_exons'][position+1]]['genomic_start'] - \
                                (self.transcriptdict['pad']*2):
                            clash_after = True
                    except KeyError:
                        print('potential undetected clash after exon {}'.format(exon_number))
                if exon_number > 1:
                    if exon_number > int(latex_dict['list_of_exons'][position-1]):
                        if ex_start < latex_dict['exons'][latex_dict['list_of_exons'][position-1]]['genomic_end'] + \
                                (self.transcriptdict['pad']*2):
                            clash_before = True
                if clash_after is True and clash_before is True:
                    self.output_list.append('BE AWARE: Flanking intron is shared with both adjacent exons')
                elif clash_after is True:
                    self.output_list.append('BE AWARE: Flanking intron is shared with the following exon')
                elif clash_before is True:
                    self.output_list.append('BE AWARE: Flanking intron is shared with the previous exon')

            sequence = exon_dict['sequence']
            characters_on_line = 0
            self.output_list.append('')
            pdfannotation_timer = 0
            for base_position in range(len(sequence)):

                # Stop each line at a specific length
                if characters_on_line % 60 == 0 and characters_on_line != 0\
                        and pdfannotation_timer == 0:
                    amino_was_printed = amino_string
                    amino_was_numbered = amino_number_string

                    if lines_on_page >= 41:
                        extra_lines = 2   # Base and numbering strings as default
                        if amino_was_numbered:
                            extra_lines += 1
                        if amino_was_printed:
                            extra_lines += 1

                        if lines_on_page + extra_lines >= 45:
                            if self.write_as_LaTex:
                                self.print_exon_end()
                            else:
                                self.output_list.append('\n\n')
                            lines_on_page = 0
                    wait_value = 0
                    amino_wait = 0
                    self.output_list.append(''.join(number_string))
                    if self.line_break_print:
                        dna_string.append('}')
                    self.output_list.append(''.join(dna_string))
                    lines_on_page += 2
                    if amino_was_printed:
                        self.output_list.append(''.join(amino_string))
                        lines_on_page += 1
                    if amino_was_numbered:
                        self.output_list.append(''.join(amino_number_string))
                        lines_on_page += 1
                    self.output_list.append('')
                    characters_on_line = 0
                    lines_on_page += 1
                    amino_string = []
                    number_string = []
                    if self.line_break_print:
                        dna_string = ['\\hl{']
                        self.line_break_print = False
                    else:
                        dna_string = []
                    amino_number_string = []
                    self.exon_spacing = False
                    self.amino_spacing = False

                char = sequence[base_position]
                if pdfannotation_timer > 0:
                    pdfannotation_timer -=1
                    dna_string.append(char)
                    if pdfannotation_timer == 0:
                        self.line_break_print = True
                elif char == '}':
                    dna_string.append(char)
                    self.line_break_print = False
                    pass

                # Deal with the insertions of PDF annotations and highlighting
                elif char == '\\':
                    dna_string.append(char)
                    subseq = sequence[base_position:]
                    match = re.search(self.pattern, subseq)
                    try:
                        pdfannotation_timer = len(match.group())-1
                    except AttributeError:
                        print(subseq)

                else:
                    dna_string.append(char)

                    if char.isupper():
                        self.exon_printed = True
                    if cds_count == 0:
                        self.amino_printing = True
                        cds_count = 1
                    if amino_acid_counter >= len(protein):
                        self.amino_printing = False
                    # Calls specific methods for character decision
                    # Simplifies local logic
                    (next_amino_string, codon_count, amino_acid_counter,
                     codon_numbered) = self.decide_amino_string_character(char, codon_count, amino_acid_counter,
                                                                          codon_numbered, protein)
                    amino_string.append(next_amino_string)
                    if next_amino_string == '*':
                        self.check_AA = False
                    if next_amino_string != ' ' and self.check_AA:
                        pos1 = char
                        check_position = base_position + 1
                        check_sequence = sequence

                        # This should only fail on the final exon; where it is not called
                        try:
                            check_next_exon = latex_dict['list_of_exons'][position+1]
                        except IndexError:
                            pass
                        if check_sequence[check_position].isupper():
                            pos2 = check_sequence[check_position]
                            check_position += 1
                        else:
                            check_sequence = latex_dict['exons'][check_next_exon]['sequence']
                            check_position = 0
                            pos2 = check_sequence[check_position]
                            while pos2.islower():
                                check_position += 1
                                pos2 = check_sequence[check_position]
                            check_position += 1
                        if check_sequence[check_position].isupper():
                            pos3 = check_sequence[check_position]
                        else:
                            check_sequence = latex_dict['exons'][check_next_exon]['sequence']
                            check_position = 0
                            pos3 = check_sequence[check_position]
                            while pos3.islower():
                                check_position += 1
                                pos3 = check_sequence[check_position]
                        
                        if pos1 == '\\' or pos1 == '}':
                            pass
                        elif pos2 == '\\' or pos2 == '}':
                            pass
                        elif pos3 == '\\' or pos3 == '}':
                            pass
                        else:
                            index = pos1+pos2+pos3
                            try:
                                if self.codon_table[index] != next_amino_string:
                                    print(
                                        "There is an error with the amino acid - codon pairing in exon {}\n"
                                        "{} - {}, Amino Acid number {}\n"
                                        "Base 3 position = {}\n"
                                        "Next few bases: {}\n".format(
                                            check_next_exon,
                                            index,
                                            next_amino_string,
                                            amino_acid_counter,
                                            check_position,
                                            check_sequence[check_position + 1:check_position + 5])
                                    )
                            except KeyError:
                                raise KeyError("The key '{}' does not have a codon entry: {}. DNA seq: {}".format(
                                    index,
                                    self.transcriptdict['genename'],
                                    dna_string
                                ))

                    (next_amino_number, amino_wait, codon_numbered,
                     amino_acid_counter) = self.decide_amino_number_string_character(amino_wait, codon_numbered,
                                                                                amino_acid_counter)
                    amino_number_string.append(next_amino_number)

                    (next_number_string, wait_value, cds_count, amino_acid_counter, post_protein_printer, intron_offset,
                     intron_in_padding, intron_out) = self.decide_number_string_character(char, wait_value, cds_count,
                                                                                          amino_acid_counter,
                                                                                          post_protein_printer,
                                                                                          intron_offset, intron_in_padding,
                                                                                          len(protein), intron_out)
                    number_string.append(next_number_string)
                    characters_on_line += 1

            # Section for incomplete lines (has not reached line-limit print)
            # Called after exon finishes printing bases
            if len(dna_string) != 0:
                if lines_on_page >= 44:
                    if self.write_as_LaTex:
                        self.print_exon_end()
                    else:
                        self.output_list.append('\n\n')
                wait_value = 0
                amino_wait = 0
                if number_string:
                    self.output_list.append(''.join(number_string))
                self.output_list.append(''.join(dna_string))
                if amino_string:
                    self.output_list.append(''.join(amino_string))
                if amino_number_string:
                    self.output_list.append(''.join(amino_number_string))
            if self.write_as_LaTex:
                self.print_exon_end()
            lines_on_page = 2
                
        for version in self.list_of_versions:
            assert isinstance(version, str)
            self.output_list.append(version)

        if self.write_as_LaTex:
            self.print_latex_footer()

    def print_exon_end(self):
        self.output_list.extend(
            [
                '\\end{alltt}',
                '\\newpage',
                '\\begin{alltt}'
                ]
        )

    def print_latex_footer(self):
        """
        A brief function to set the final lines of the document if the output
        is to be in a LaTex parse-able format
        """
        self.output_list.extend(
            [
                '\\end{alltt}',
                '\\end{document}'
                ]
        )

    def line_printer(self, string):
        """
        :param string: next string to be written to output list
        :return: none

        Generic print method to handle all list output in one location
        """
        self.output_list.append(''.join(string))

    @staticmethod
    def decide_amino_number_string_character(amino_wait, codon_numbered, amino_acid_counter):
        output = ''
        if amino_wait != 0:
            amino_wait -= 1
        elif amino_wait == 0:
            if amino_acid_counter % 10 == 1:
                if not codon_numbered:
                    output = '|' + str(amino_acid_counter)
                    amino_wait = len(str(amino_acid_counter))
                    codon_numbered = True
                elif codon_numbered:
                    output = ' '
            elif amino_acid_counter % 10 != 1:
                output = ' '
        return output, amino_wait, codon_numbered, amino_acid_counter

    def run(self):
        self.print_latex()
        return self.output_list, self.nm
