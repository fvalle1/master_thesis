# Network analysis on TCGA data


# Files
## TCGA
### TCGA get manifest
[TCGA_GetManifest](TCGA_GetManifest.ipynb) creates a file called *manifest.txt* that is useful to download data using `gdc-client download -m manifest.txt` command

### TCGA API
[TCGA_API](TCGA_API.ipynb) is useful to retrieve informations about a single file (sample)

##Table
### Table creation
[Table_Creation](Table_Creation.ipynb) reads a folder with the data downloaded with **gdc-client** and creates a mainTable.csv dataset

### Table mining
[Table_Mining](Table_Mining.ipynb) is useful to extract means and vars per tissues. It stores a *meanVariances.csv* file useful in next step

### Table Analyzer
[Table_Analyzer](Table_Analyzer.ipynb) is useful to plot expression per tissues histograms and FPKM means and distributions. Here it mean versus occurrence is created.

### Table by gene type
[Table_ByGeneType](Table_ByGeneType.ipynb) is useful to do analysis separating protein coding from non coding excetera

### Table by tissue
[Table_ByTissue](Table_ByTissue.ipynb) perform analysis looking at different tissues

### Table protein coding
[Table_ProteinCoding](Table_ProteinCoding.ipynb) isolate a *mainTable.csv* file with only protein coding genes. Creates a *genes_pc.txt* file with all ENSG-id that are related to a protein-coding gene

#### Table microRNA
[Table_microRNA](Table_microRNA.ipynb) as above but with microRNA

### Table Analyzer
[Table_tfidf](Table_tfidf.ipynb) perform tf-idf analysis


## hSBM
### topics
[hSBM_topics](hSBM_topics.ipynb) analyse topic composition

### clusters
[hSBM_clusters](hSBM_clusters.ipynb) bench the performance in clustering samples

### topic_dist
[hSBM_topic-dist](hSBM_topic-dist.ipynb) analyse topic distribution inside samples



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
0 ---> read mainTable.csv
1 ---> read and extimate correlation mainTable.csv
2 ---> extimate means and variances
3 ---> GenerateNullData
4 ---> read nullTable.csv
5 ---> nullTable.csv read and extimate correlation
6 ---> nullTable.csv extimate means and variances
7 ---> read and make bipartite graph
 0.811598s wall, 0.160000s user + 0.010000s system = 0.170000s CPU (20.9%)
```

Options with **read** are able to create
* *A.dat* file with **abundances**
* *O.dat* file with **occurrences**
* *heaps.dat* file with **sizes** and **vocabulary sizes**

Options with **extimate correlations** create a *correlations.dat* file with **H(X:Y)** for each couple of words.

**Generate null data** option create a *nullTable.csv* file with null model generated data.

**Make bipartite graph** option creates a *graph.xml.gz* file which is a great input *hieratical Stochastic Block Model*

## Zipf
Use [Zipf.ipynb](Zipf.ipynb) to plot Zipf and U (occurrence distribution)

## Heaps
Use [Heaps.ipynb](Heaps.ipynb) to plot Heaps distibution
