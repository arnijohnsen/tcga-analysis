# Introduction
This repository contains R scripts which analyse data from 
The Cancer Genome Atlas:

- Breast invasive carcinoma (brca)
- Ovarian serous cystadenocarcinoma (ov)
- Pancreatic adenocarcinoma (paad)
- Prostate adenocarcinoma (prad)

Data is available for download at https://tcga-data.nci.nih.gov/tcga/

Author: Arni Johnsen, arni.johnsen@gmail.com

# Installation 

Installation: 

    git clone https://github.com/arnijohnsen/tcga-analysis.git
    mkdir raw-data
    mkdir parsed-data

The Rscripts rely on the following packages:

    data.table
    limma

# Examples

See EXAMPLE.md (none available currently)

# Structure

## Directory structure

Each cancer type (e.g. brca) there are several data directories:

    .
    ├── parsed-data
    │   └── brca
    │       ├── cnv
    │       ├── expr
    │       ├── info
    │       ├── meth
    │       ├── mirn
    │       └── muta
    └── raw-data
       ├── annotation-files
       └── brca
           ├── clin
           ├── cnv
           ├── expr
           ├── meth
           ├── mirn
           └── muta

tar files from TCGA should be untarred to one of

- `rawdata/brca/cnv` (copy number variation)
- `rawdata/brca/expr` (RNA sequencing)
- `rawdata/brca/meth` (DNA methylation)
- `rawdata/brca/mirn` (miRNA sequencing)
- `rawdata/brca/muta` (somatic mutations)

File reading scripts will save .Rdata files to subdirectories of `parsed-data`,
which contain data.tables with data from files in `raw-data`.

## Data table structure

Data contained in `pased-data/xxxx/yyyy` has a consistant
nomenclature and format:

- Data files are name `xxxx-yyyy-zzzzzz.Rdata`, where
    - `xxxx` is the cancer type (e.g. brca)
    - `yyyy` is the type of data (e.g. expr)
    - `zzzzzz` is either `cancer` or `normal`
- Each data file contains one `data.frame`, named `xxxx.yyyy.zzzzzz`, following the same nomenclature
- The first column(s) of each data.table is the name of probe or gene
- Each column of a data frame represents one sample and is named by its TCGA barcode

TODO: Update this description, as it is different for mutation and copy number data
