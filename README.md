# ExAC_nonTCGA_explore_v6
A shiny app that opens filtered ExAC non-TCGA data on a web browser for exploring

Can be accessed with the following link:
(https://medgen.shinyapps.io/exac_nontcga_explore_v6/)

Provides access to a filtered ver0.3.1 ExAC dataset containing only Loss of function, inframe indels, and predicted deleterious and damaging missense variants (by SIFT and PolyPhen respectively).

# Additional filters
Allows for additional filtering on predicted variant consequence:
## Variant impact:
### Protein affecting (default)
- Loss of function
- Inframe indels
- predicted deleterious and damaging missense variants (by SIFT and PolyPhen respectively)

### Loss of function
- Stop gained
- Stop lost
- Start lost
- Frameshift variant
- Splice donor and splice acceptor variants

## Variant allele frequency
Allows for additional filtering on variant rarity:
### All (default)
### < 0.05
Includes all variants with an allele frequency of < 0.05 within the entire ExAC non-TCGA set based on the AF column
### < 0.01
Includes all variants with an allele frequency of < 0.01 within the entire ExAC non-TCGA set based on the AF column

## Gene
Allows the input of incomplete strings and will return all genes within the data that contain that string. Has the option to match the gene name exactly from the search. 

# Data output
Large table by default contains all protein affecting variants and is updated with filters. Hitting the download button downloads all data that match filters. It is currently not possible to download data without first selecting a gene of interest.

# Summary Metrics
Summary metrics are generated after a gene name is input into the set. Counts of occurences (sum of ACs) are summarised for different impact and AF variants for each ExAC population. If searched gene matches several genes then all genes will be summarised unless otherwise specified with the 'Match gene name exactly' button.
