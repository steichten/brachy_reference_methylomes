#Analysis scripts for _Brachypodium distachyon_ reference methylomes
---

##This provides all of the scripts required to reproduce the results of:

Eichten SR, Stuart T, Srivastava A, Lister R, Borevitz JO. DNA methylation
profiles of diverse Brachypodium distachyon aligns with underlying genetic
diversity. Genome Res. 2016 Sep 9. pii: gr.205468.116. [Epub ahead of print]
PubMed PMID: 27613611.
---

#Software requirements

In order to run all of these you will need:

- A fair bit of working HD space (think 400GB+)
- SRA tools (https://github.com/ncbi/sra-tools/wiki/Downloads)
- bedtools (https://github.com/arq5x/bedtools2)
- bismark v0.13.0 (http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- trimgalore (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- bowtie1 (http://bowtie-bio.sourceforge.net/index.shtml)
- R (https://www.r-project.org/)

In addition, there are many R packages required:

- ggplot2
- plyr
- dplyr
- fields

