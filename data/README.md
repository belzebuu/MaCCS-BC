The format must be such that the amount of work required to put the data
in the required form is minimized with respect to what available from:
http://cbio.mskcc.org/cancergenomics/pancan_tcga/


- a file for the interaction network (ie, undirected graph) containing a
  list of edges one per row:

  For example, content of file BRCA.txt:
  1433B 1433E 
  1433B 1A1L1
  ...

  It indicates that 1433B, 1433E, 1A1L1 are identifiers of vertices of the
  interaction graph and that (1433B, 1433E) and (1433B, 1A1L1) are two
  edges. Edges can be written just once or twice, the parser should
  take care of recognizing that something has already been expressed: 
   
  Eg: 
   1433B 1433E 
   1433E 1433B 
   ...
  are the same.


  Alternatively, one can give the list of adjacent vertices:
  
  1433B 1433E 1A1L1

  In this case the format can be redundant, that is, 1433B will appear
  also when the adjacent vertices of 1433E are listed. Each list ends
  with a new line.



- coverage network. Here I would try to stay as close as possible to the
  format of http://cbio.mskcc.org/cancergenomics/pancan_tcga/

  The file contains a collection of adjacency lists, one per row. Each
  adjacency list represents the patients that had a specific gene
  mutated.
  
  The first term of each row indicates the type of cancer, the second
  the name of the gene and the following terms indicate the patients
  that had the gene mutated.
 
  [Please indicate here if the role of genes and patients is exchanged]
  
  For example, the file laml_pancancer.mm contains:

  laml    TCGA-AB-2802    ANKRD30A        C20orf24        C6orf10 CYP21A2 DNMT3A  ESYT1   IDH1    KIAA2022        KRT18P40
        LOC100134184    LOC728135       LOC728211       PTPN11  TBX15   TCHHL1
  laml    TCGA-AB-2803    ABCC1   ASMTL   CACNA1S CC2D2A  ENSG00000222666 FNIP2   FSD2    HMX1    KDM2B   LMOD1   LNX1    LOC100131868    LOC100132516    LOC100133234    LOC147804       LOC643342       MPPE1   MT-CO1  PEX13   POU5F1  RNPEPL1 RPL32   RYR3    TRIM48  TUBA1A  WNK1
  ...
  
  
  which indicates that in cancer laml the gene TCGA-AB-2802 was found
  mutated in the patients listed.
  
  [do we need to have this file containing indication of the different type of cancer? I th ink not, a user could simply keep different cancers in different files and concatenate them in a single file when needed. Hence, I would remove the first identifier for the cancer]
 
 
  
  
- weight file. This file is optional. It contains a list of weights
  indicating the importance of genes and/or edges. Weights must be from
  the interval [0,1].

  Each row of the file contains a gene or an edge plus the corresponding
  weight.

  For example:
  A1BG    0.0
  A1CF    0.270809121012
  A2BP1   0.00541650013642
  A2M     0.6426993429
  A1BG KIAA2022 0.876
  A1BG ASMTL 0.876
  ...
 
  it indicates that all connections from gene A1BG to patients have
  weight 0.0 except for the connections to patients KIAA2022 and ASMTL
  that have the weights indicated. 

  Genes not listed are assigned weight 0.



The parsers must handle both space or tab separated lines.
