from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

samples = pd.read_csv("files.dat", index_col=[0])

def get_generic_tissue_from_specific(tissue, samples=samples):
    return samples[samples['secondary_site']==tissue]['primary_site'].values[0]

def get_specific_mapping_to(tissue, samples=samples):
    return samples[samples['primary_site']==tissue]['secondary_site'].unique()

def get_gtex_tissue(sample):
    for fullsample in samples.index.values:
        if sample in fullsample:
            return samples.loc[fullsample,:]
        
def makePie(df_clusters, level, c, whatToLookFor = ['primary_site','secondary_site']):
    c=int(c)
    level=int(level)
    fig = plt.figure(figsize=(60,15))
    ax = fig.subplots(1, len(whatToLookFor))
    for i,lookFor in enumerate(whatToLookFor):
        datatotestarr = []
        for sample in df_clusters['Cluster %d'%c].dropna():
            try:
                datatotestarr.append(get_gtex_tissue(sample)[lookFor])
            except:
                print("error with %s"%sample)
        utype, counts = np.unique(datatotestarr, return_counts=True)
        total = len(datatotestarr)
        try:
            labels = ['\n'.join(wrap(str(l), 20)) for l in utype]
        except:
            labels = utype
        ax[i].set_title(lookFor, fontsize=44)
        patches, texts, autotexts = ax[i].pie(counts,
                                              labels=labels,
                                              autopct=lambda p: '#:%.0f'%(p * total / 100),
                                              textprops={'fontsize':30, 'color':'white', 'wrap':True})
        for t in texts:
            t.set_fontsize(24)
            t.set_wrap(True)
            t.set_color('black')
    fig.savefig("cluster_pie_level_%d_cluster_%d.png"%(level, c))
    
def get_cluster_given_l(l, directory):
    df_clusters = pd.read_csv("%s/topsbm/topsbm_level_%d_clusters.csv"%(directory, l), header=[0])
    cluster={}
    for i,c in enumerate(df_clusters.columns):
        cluster[i]=df_clusters[c].dropna().values
    return cluster
    
def define_labels(cluster, label='primary_site', verbose=False):
    true_labels = []
    predicted_labels = []
    for c in cluster:
        if verbose:
            print(c)
        for sample in cluster[c]:
            #true_labels.append(getFile(sample)['primary_site'].values[0])
            try:
                true_labels.append(get_gtex_tissue(sample)[label])
                predicted_labels.append(c)
            except:
                print("error in %s"%sample)
    _, true_labels = np.unique(true_labels,return_inverse=True)
    return true_labels, predicted_labels


def get_fraction_sites(cluster, df_files, label='primary_site', normalise=False):
    fraction_sites = {}
    c_fraction_site = {}
    for site in np.unique(df_files[label].values):
        fraction_sites[site] = []
        c_fraction_site[site] = 0

    for i,c in enumerate(cluster):
        for sample in cluster[i]:
            try:
                for fullsample in df_files.index.values:
                    if sample in fullsample:
                        foundsample=df_files.loc[fullsample,:]
                c_fraction_site[foundsample[label]]+=1
            except:
                print("error in %s"%sample)
        for site in fraction_sites:
            if normalise:
                norm = float(len(cluster[i]))
            else:
                norm = 1
            fraction_sites[site].append(c_fraction_site[site]/norm)
            c_fraction_site[site]=0
    return fraction_sites


def get_clustersinfo(cluster, fraction_sites):
    clustersinfo = {
    "maximum": [],
    "sizes": [],
    "nclasses":[]
    }
    for icluster in cluster:
        maximum = 0
        size = 0
        nclass = 0
        site_maximum = ''
        cumulative = 0
        for site, data in fraction_sites.items():
            cdata = data[icluster]
            cumulative += cdata
            if cdata > maximum:
                maximum = cdata
                site_maximum = site
            if cdata > 0:
                nclass+=1
            size+=cdata
        clustersinfo['maximum'].append([float(maximum)/cumulative,site_maximum])
        clustersinfo['sizes'].append(size)
        clustersinfo['nclasses'].append(nclass)
        maximum=0
        cumulative=0
        size=0
        nclass=0
        site_maximum=''
    return clustersinfo