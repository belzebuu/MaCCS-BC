# Data Format


- A file for the interaction network (ie, undirected graph) containing a
  list of edges one per row:

  For example, content of file BRCA.txt:
  1433B 1433E 
  1433B 1A1L1
  ...

  It indicates that 1433B, 1433E, 1A1L1 are identifiers of vertices of the
  interaction graph and that (1433B, 1433E) and (1433B, 1A1L1) are two
  edges. Edges can be written just once or twice, the parser 
  takes care of recognizing that something has already been expressed: 
   
  Eg: 
   1433B 1433E 
   1433E 1433B 
   ...
  are the same.



- A file for the coverage network. The format is similar to the one at
  http://cbio.mskcc.org/cancergenomics/pancan_tcga/

  The file contains a collection of adjacency lists, one per row. Each
  adjacency list represents the patients that had a specific gene
  mutated.
  
  Each row reports the genes mutated in the same type of cancer of a
  patient. The first term of each row is an identifier for the patient
  and the terms that follow indicate the genes that were mutated in the
  cancer.
 
  For example, the file laml_pancancer.mm for the type of cancer laml contains:

  TCGA-AB-2802    ANKRD30A        C20orf24        C6orf10 CYP21A2 DNMT3A  ESYT1   IDH1    KIAA2022        KRT18P40
        LOC100134184    LOC728135       LOC728211       PTPN11  TBX15   TCHHL1
  TCGA-AB-2803    ABCC1   ASMTL   CACNA1S CC2D2A  ENSG00000222666 FNIP2   FSD2    HMX1    KDM2B   LMOD1   LNX1    LOC100131868    LOC100132516    LOC100133234    LOC147804       LOC643342       MPPE1   MT-CO1  PEX13   POU5F1  RNPEPL1 RPL32   RYR3    TRIM48  TUBA1A  WNK1
  ...
  
  
  which indicates that the patient TCGA-AB-2802 had the genes listed
  mutated in the cancer of type laml.
    
  
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

  Genes and edges not listed are assigned weight 0.



The parser interprets both space or tab separated lines.
