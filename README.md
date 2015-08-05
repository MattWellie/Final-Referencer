##License 
At some point I'll work out what kind of license this should have...
It's open source anyway...

##Requirements and testing

pdflatex/texlive/miktex (For LaTex production, not required for text document (see write\_as\_latex in Xml_GUI)
-   This requires the packages Soul & Color (Default LaTex packages), Pdfcomment (possibly a default package, https://www.ctan.org/pkg/pdfcomment?lang=en), and alltt (available by default)

tkinter (python graphics package)

BioPython 1.65 (GenBank file parsing, not required for LRG file references)

Python 2.7

##What it does

- This program takes a file as input and creates a typeset reference sequence.
This sequence is designed to be used during sequence checking, and includes 
HGVS nomenclature base and amino acid numbering

- This program ensures platform independence through use of OS-ambivalent file paths
and naming conventions. 

- Separation into discrete paths allows files in valid GenBank or LRG formats to be used as input
without making any changes to the program code.

- GenBank format files sourced from Ensembl are also usable, with a set of conditions
    present to handle the format differences between NCBI and Ensembl

##Who its for

This is intended for use by any genetic scientist who may require a hard-copy of a genetic sequence 
for checking against when signing off reports. This is intended to replace the existing procedure
of manually creating reference sequence copies by using a uniform typesetting style and format for
all inputs.

##How to operate it
- At command line: *python XML_GUI.py* without any additional arguments. This can be done using *python X -> Tab*

- This command will bring up the main interface, offering two entry boxes. 

- The top entry box can either be edited directly or by using the 'Browse...' button. This 
will show the local directory and allow file selection directly. This is where to insert the
filename you wish to convert. The default contents of this box can be set in XML_GUI.py.

- There is a 'Help' button on the ribbon which will print a brief guide statement
to the command line. This can be edited in the XML_gui.py file

##How it works

- The program has been broken up into several different components;
    - XML_gui.py to show the user interface
    - LRG/GBK_Parser.py to read the input file into a dictionary
    - optional call to primer module to annotate primers in final output
    - reader.py to read the dictionary into a list output format
    - writer.py to read the list into an actual file
    - latex_writer.py to write the reader output into a external file 
    - The XML_GUI.py module then calls a pdflatex command to typeset the file

- For GB and LRG files with multiple transcripts the program has separate ways of dealing with contents
    - For .gb files from NCBI, the program will only use CDS and mRNA features which have a gene 
        annotation matching the gene name attached to each Exon. For files sourced from NCBI, the 
        file may contain multiple genes and isoforms which span the region of the main gene, but 
        only the main gene name will feature in exon annotations
    - For .gb files from Ensembl, manual hacking is required. After downloading a copy of the Ensembl 
        output in .gb format (Export Data -> Output:GB -> uncheck all but 'Gene Information'), 
        delete all CDS and mRNA features which do not correspond to the accession of choice (usually 
        all but the first pair), and try-catch blocks will handle the rest of the processing. This can
        be used as standard GB input once the changes are made. Less manual solution being investigated.
    - LRG files do not contain details of other genes spanning the region, so each of the separate <transcript>
        blocks is handled independently, along with the corresponding sets of exon coordinates. This offers the 
        same content as the GB files, though the format is clearer for parsing.
- In all cases, if multiple valid transcripts exist, a separate file is printed for each. These are currently numbered
    sequentially (1, 2...) and will contain the specific transcript details within the file contents

##Planned updates

- No further updates planned

##Issue Tracker

- Kludgey fixes in place in Gbk_Parser.py to allow for use of Ensembl transcripts in gb format, which do not
    contain same feature annotations as NCBI GenBank files. This may cause as-yet undetected problems, check 
    output

###Tinkering:

To have the amino acid appear below different bases within the codon:
* set the variable codon_count in print_latex() to:
    - 3 to print with the first base of the codon
    - 2 to print in the middle
    - 1 to print below the final base

The protein sequence is grabbed as a single uninterrupted sequence, so none of these 
approaches will affect the number of AAs printed, only the positions. However, this 
may change which exon the AA is printed in (if a codon is across exons)

This program makes use of the ArgParse Python module (25/03/2015). The optional arguments listed can be 
used to alter the function of the program:
* --trim : prevent intronic flanking sequence being trimmed to prevent overlapping sequence
* --clash : prevent messages being printed under exon headers to indicate sequence overlap
* --text : print output files as text only (rather than being processed with LaTex)


─────────▄──────────────▄<br>
────────▌▒█───────────▄▀▒▌<br>
────────▌▒▒▀▄───────▄▀▒▒▒▐<br>
───────▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐<br>
─────▄▄▀▒▒▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐<br>
───▄▀▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀██▀▒▌<br>
──▐▒▒▒▄▄▄▒▒▒▒▒▒▒▒▒▒▒▒▒▀▄▒▒▌<br>
──▌▒▒▐▄█▀▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐<br>
─▐▒▒▒▒▒▒▒▒▒▒▒▌██▀▒▒▒▒▒▒▒▒▀▄▌<br>
─▌▒▀▄██▄▒▒▒▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌<br>
─▌▀▐▄█▄█▌▄▒▀▒▒▒▒▒▒░░░░░░▒▒▒▐<br>
▐▒▀▐▀▐▀▒▒▄▄▒▄▒▒▒▒▒░░░░░░▒▒▒▒▌<br>
▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒░░░░░░▒▒▒▐<br>
─▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌<br>
─▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▐<br>
──▀▄▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▄▒▒▒▒▌<br>
────▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀<br>
───▐▀▒▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀<br>
 ──▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀▀<br>
