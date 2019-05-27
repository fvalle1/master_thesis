import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def plot_cluster_composition(fraction_sites, directory, level, normalise=False, label='primary_site', shuffled=False):
    sns.set(font_scale=0.8)
    df_clusters = pd.read_csv("%s/topsbm/topsbm_level_%d_clusters.csv"%(directory, level), header=[0])
    x = np.arange(1,1+len(df_clusters.columns))
    bottom = np.zeros(len(x))
    ymax = 0
    fig=plt.figure(figsize=(15,8))
    for site, data in fraction_sites.items():
        if np.max(data) == 0:
            continue
        plt.bar(x,data,label=site, bottom=bottom)
        bottom=bottom+data
    plt.xlabel("cluster", fontsize=16)
    if normalise:
        plt.ylabel("fraction of %s"%label, fontsize=16)
    else:
        plt.ylabel("# of %s"%label, fontsize=16)
    plt.title("%s distribution across clusters"%label, fontsize=16)
    plt.legend(ncol=3)
    plt.xticks(x)
    plt.show()
    fig.savefig("%s/%s%sclustercomposition_l%d_%s.png"%(directory, "shuffled" if shuffled else '', "fraction_" if normalise else '', int(level),label))

def get_cluster_given_l(l, directory):
    df_clusters = pd.read_csv("%s/topsbm/topsbm_level_%d_clusters.csv"%(directory, l), header=[0])
    cluster={}
    for i,c in enumerate(df_clusters.columns):
        cluster[i]=df_clusters[c].dropna().values
    return cluster

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

def plot_maximum(clustersinfo, cluster, label, level, directory,clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(15,6))
    real = np.array(clustersinfo['maximum'])[:,0].astype(float)
    plt.plot(np.sort(real), marker='x')
    shuffled=False
    if clustersinfo_shuffle is not None:
        shuffled = np.array(clustersinfo_shuffle['maximum'])[:,0].astype(float)
        plt.plot(np.sort(shuffled), marker='x')
        plt.plot(np.sort(real-shuffled), marker='x')
        shuffled=True
    plt.plot(np.arange(len(cluster)),[0.8 for i in range(len(cluster))], visible=True)
    plt.xlabel("cluster", fontsize=16)
    plt.ylabel("maximum fraction\nwith same %s"%label, fontsize=16)
    plt.ylim((0,1.1))
    plt.show()
    fig.savefig("%s/%scluster_maximum_l%d_%s.png"%(directory,"shuffled" if shuffled else '',level,label))



def plot_maximum_size(clustersinfo, label, level, directory,clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(15,6))
    x = np.array(clustersinfo['sizes']).astype(int)
    y = np.array(clustersinfo['maximum'])[:,0].astype(float)
    plt.scatter(x,y, lw=10, label='clusters')
    plt.xlim(0,np.max(x)+np.max(x)/10)
    plt.plot([10 for _ in range(10)],np.linspace(0,1,10), ls='--', label=10)
    plt.plot(np.linspace(0.5,x.max()),1./np.linspace(0.5,x.max()), label='uniform')
    shuffled=False
    if clustersinfo_shuffle is not None:
        shuffled=True
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:,0].astype(float)
        plt.scatter(x_shuffle,y_shuffle, lw=10, label='clusters')
        plt.xlim(0,np.max(x_shuffle)+np.max(x_shuffle)/10)
    plt.xlabel("cluster size", fontsize=16)
    plt.ylabel("maximum fraction\nwith same %s"%label, fontsize=16)
    plt.ylim((0,1.1))
    plt.legend(loc='lower right', fontsize=18)
    plt.show()
    fig.savefig("%s/%sclusterhomosize_l%d_%s.png"%(directory, "shuffled" if shuffled else '', level,label))

def plot_maximum_label(clustersinfo,label,level, directory,clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(10,6))
    x = np.array(clustersinfo['nclasses']).astype(int)
    y = np.array(clustersinfo['maximum'])[:,0].astype(float)
    shuffled=False
    plt.scatter(x,y, lw=10, alpha=0.5, label='clusters')
    plt.plot(np.arange(1,np.max(x)+2),1./np.arange(1,np.max(x)+2),ls='--',c='cyan', label='uniform')
    plt.xlim(0.95,np.max(x)+0.5)
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:,0].astype(float)
        plt.scatter(x_shuffle,y_shuffle, lw=10, alpha=0.5, label='clusters shuffled')
        plt.plot(np.arange(1,np.max(x_shuffle)+2),1./np.arange(1,np.max(x_shuffle)+2),ls='--',c='cyan', label='')
        shuffled=True
        plt.xlim(0.95,np.max(x_shuffle)+0.5)
    plt.xlabel("# labels in cluster", fontsize=16)
    plt.ylabel("maximum fraction\nwith same %s"%label, fontsize=16)
    plt.ylim((0,1.1))
    plt.legend(loc='lower right', fontsize=18)
    plt.show()
    fig.savefig("%s/%scluster_homon_l%d_%s.png"%(directory, "shuffled" if shuffled else '', level,label))

def plot_labels_size(clustersinfo, label, level,directory, clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(10,6))
    x = np.array(clustersinfo['sizes']).astype(float)
    y = np.array(clustersinfo['nclasses']).astype(int)
    plt.xlim(x.min()-10,x.max()+5)
    plt.ylim(y.min()-2,y.max()+5)
    shuffled=False
    plt.scatter(x,y, lw=10, alpha=0.5, label='clusters')
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(float)
        y_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        plt.scatter(x_shuffle,y_shuffle, lw=10, alpha=0.5, label='clusters shuffled')
        plt.xlim(x.min()-10,x_shuffle.max()+5)
        plt.ylim(y.min()-2,y_shuffle.max()+8)
        shuffled=True
    plt.xlabel("cluster size", fontsize=16)
    plt.ylabel("# labels in cluster", fontsize=16)
    plt.legend(loc='upper right', fontsize=18)
    plt.show()
    fig.savefig("%s/%scluster_shuffle_label_size_l%d_%s.png"%(directory, "shuffled" if shuffled else '', level,label))

def make_heatmap(fraction_sites, directory, label, level, shuffled=False, normalise=False):
    sns.set(font_scale=2)
    found_classes = []
    for site, data in fraction_sites.items():
        if np.max(data) == 0:
            continue
        found_classes.append(site)
    for arr in fraction_sites.values():
        x = len(arr)
        break
    x = np.arange(1,1+x)
    fig = plt.figure(figsize=(30,10))
    heatmap = fig.subplots(1)
    heatmap = sns.heatmap(pd.DataFrame(data=fraction_sites).loc[:,found_classes].transpose(), vmin=0, cmap= "RdYlBu_r", xticklabels=x)
    fig.savefig("%s/%sheatmap_cluster%s_l%d_%s.png"%(directory,"shuffled" if shuffled else '',"fraction_" if normalise else '', int(level),label))

def get_file(sample, df_file):
    for fullsample in df_file.index.values:
        if sample in fullsample:
            return df_file.loc[fullsample,:]

def define_labels(cluster, df_files, label='primary_site', verbose=False):
    true_labels = []
    predicted_labels = []
    for c in cluster:
        if verbose:
            print(c)
        for sample in cluster[c]:
            try:
                true_labels.append(get_file(sample, df_files)[label])
                predicted_labels.append(c)
            except:
                print("error in %s"%sample)
    _, true_labels = np.unique(true_labels,return_inverse=True)
    return true_labels, predicted_labels
