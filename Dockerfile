FROM cs-prod-tools-artifactory-01.gel.zone:5004/bertha/bertha-compute

RUN pip install biopython==1.65 && sudo apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
