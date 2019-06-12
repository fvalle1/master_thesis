# Network analysis on TCGA data

## Pipeline
## Run statistics
```python
import tacos
tacos.statistics()
tacos.tacos.shuffleLabels()
```

### make bipartite network
```bash
python addo.py
```
```python
import tacos
tacos.makegraph()
```

## Run stochastic block model
```
docker run -it -p8888:8888 -v ~/home/cloudadm/drive/:/home/user/ fvalle01/hsbm bash
jupyter notebook --ip0.0.0.0 --allow-root
```

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
Use [Tool for Analyse COmponents Systems](tacos)

## Zipf
Use [Zipf.ipynb](Zipf.ipynb) to plot Zipf and U (occurrence distribution)

## Heaps
Use [Heaps.ipynb](Heaps.ipynb) to plot Heaps distibution

# Topic Modelling

Topic model can be run from Docker repository: https://hub.docker.com/r/fvalle01/hsbm
```bash
docker run -p8888:8888 -it fvalle01/hsbm bash

jupyter notebook --ip 0.0.0.0 --allow-root
```
