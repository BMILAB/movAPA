# movAPA

Modeling and visualization of dynamics of alternative polyadenylation across biological samples

## About

Alternative polyadenylation (APA) has been widely recognized as a widespread mechanism modulated dynamically. Studies based on 3′ end sequencing and/or RNA-seq have profiled poly(A) sites in various species with diverse pipelines, yet no unified and easy-to-use toolkit is available for comprehensive APA analyses. We developed an R package called __movAPA__ for modeling and visualization of dynamics of alternative polyadenylation across biological samples. movAPA incorporates rich functions for preprocessing, annotation, and statistical analyses of poly(A) sites, identification of poly(A) signals, profiling of APA dynamics, and visualization. Particularly, seven metrics are provided for measuring the tissue-specificity or usages of APA sites across samples. Three methods are used for identifying 3′ UTR shortening/lengthening events between conditions. APA site switching involving non-3′ UTR polyadenylation can also be explored. Using poly(A) site data from rice and mouse sperm cells, we demonstrated the high scalability and flexibility of movAPA in profiling APA dynamics across tissues and single cells.

The movAPA package consists of seven main modules (Figure 1).
(1) Poly(A) sites of biological samples obtained from 3’ seq or RNA-seq are stored in the PACdataset object and further preprocessed for the removal of internal priming artifacts and normalization.
(2) The genome annotation file in GFF3/GTF format is processed, then poly(A) sites are annotated with rich information such as gene id, gene type, and genomic regions. 
(3) Statistical analyses can be conducted to profile the global landscape of poly(A) site distributions and count the overlap with other poly(A) site datasets. 
(4) Sequences surrounding poly(A) sites can be extracted and poly(A) signals of specific regions can be identified. 
(5) Three metrics can be adopted for the quantification of the usage of each poly(A) site across samples and four metrics are used for the quantification of dynamic APA site usage of a gene. (6) APA dynamics across biological samples can be profiled, including the detection of differentially expressed poly(A) sites and genes, 3′ UTR lengthening/shortening events, and canonical or non-canonical APA site switching events. (7) Rich functions are provided for the visualization of poly(A) site distributions and dynamic APA site usages across selected biological samples or single cells. 

## Getting started
### Installation


