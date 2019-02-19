# Network analysis on TCGA data


# Files
## TCGA GetManifest
[TCGA_GetManifest](TCGA_GetManifest.ipynb) creates a file called *manifest.txt* that is useful to download data using `gdc-client download -m manifest.txt` command

## Data Table
[DataTable](DataTable.ipynb) reads a folder with the data downloaded with *gdc-client* and creates a mainTable.csv dataset

## Table Analyzer
[TableAnalyzer](TableAnalyzer.ipynb) is useful to extract only protein-coding rows and to plot expression per tissues histograms

## Protein coding
[ProteinCoding](Protein Coding.ipynb) creates a *genes_pc.txt* file with all ENSG-id that are related to a protein-coding gene

## TCGA API
[TCGA_API](TCGA_API.ipynb) is useful to retrieve informations about a single file (sample)

# Analysis
## data mining
In data_mining folder run:
```
mkdir master && cd master
cmake ..
make data_mining
```

When run `./data_ming` an *help* is printed about

```bash
Running TCGA
threads: 2
Please write some options
0 ---> read and extimate correlation mainTable.csv
1 ---> read mainTable.csv
2 ---> GenerateNullData
3 ---> read nullTable.csv
4 ---> read and extimate correlation nullTable.csv
5 ---> read and makeCorpus
 0.811598s wall, 0.160000s user + 0.010000s system = 0.170000s CPU (20.9%)
```

Options with **read** are able to create
* *A.dat* file with **abundances**
* *O.dat* file with **occurrences**
* *heaps.dat* file with **sizes** and **vocabulary sizes**

Options with **extimate correlations** create a *correlations.dat* file with **H(X:Y)** for each couple of words.

**Generate null data** option create a *nullTable.csv* file with null model generated data.

**Make Corpus** option creates a *corpus.txt* file which is a great input *hieratical Stochastic Block Model*

## Zipf
Use [Zipf.ipynb](Zipf.ipynb) to plot Zipf and U (occurrence distribution)

## Heaps
Use [Heaps.ipynb](Heaps.ipynb) to plot Heaps distibution
