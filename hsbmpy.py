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
    fig.savefig("%s/%s%sclustercomposition_l%d_%s.pdf"%(directory, "shuffled" if shuffled else '', "fraction_" if normalise else '', int(level),label))

def get_cluster_given_l(l, directory):
    df_clusters = pd.read_csv("%s/topsbm/topsbm_level_%d_clusters.csv"%(directory, l), header=[0])
    cluster={}
    for i,c in enumerate(df_clusters.columns):
        cluster[i]=df_clusters[c].dropna().values
    return cluster

def get_fraction_sites(cluster, df_files, label='primary_site', normalise=False):
    fraction_sites = {}
    c_fraction_site = {}
    for site in df_files[label].dropna().unique():
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
    df = pd.DataFrame(data=fraction_sites)
    ##put first columns that have high values in average
    avgs=df.apply(lambda x: np.average(x.to_numpy()[x.to_numpy().nonzero()[0]]),axis=0)
    df=df.transpose()
    df.insert(0,'avg',avgs)
    df = df.sort_values(by=['avg'], axis=0, ascending=False).drop('avg',axis=1).transpose()
    df = df.sort_values(by=[tissue for tissue in df.columns], axis=0, ascending=False)
    return df.to_dict(orient='list')

def get_clustersinfo(cluster, fraction_sites):
    clustersinfo = {
    "maximum": [],
    "homogeneity":[],
    "sizes": [],
    "nclasses":[]
    }
    for icluster in cluster:
        maximum = 0
        homo = 0
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
                #using fraction_items normalised
                if cdata<=1:
                    homo-=cdata*np.log(cdata)
            size+=cdata
        if cumulative > 0:
            clustersinfo['maximum'].append([float(maximum)/cumulative,site_maximum])
        else:
            clustersinfo['maximum'].append([0,site_maximum])
        clustersinfo['sizes'].append(size)
        clustersinfo['nclasses'].append(nclass)
        clustersinfo['homogeneity'].append(1-homo)
        homo=0
        maximum=0
        cumulative=0
        size=0
        nclass=0
        site_maximum=''
    return clustersinfo

def plot_maximum(clustersinfo, cluster, label, level, directory,clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(15,6))
    ax = fig.subplots(1,2)
    bins = 10 
    real = np.array(clustersinfo['maximum'])[:,0].astype(float)
    ax[0].plot(np.sort(real), marker='o', ms=25, ls='')
    ax[1].hist(np.sort(real), histtype='step', bins=bins, lw=4, density=True, range=(0.05,1.05))
    shuffled=False
    if clustersinfo_shuffle is not None:
        shuffled = np.array(clustersinfo_shuffle['maximum'])[:,0].astype(float)
        ax[0].plot(np.sort(shuffled), marker='o', ls='', ms=25)
        ax[1].hist(np.sort(shuffled), histtype='step', bins=bins, lw=4, density=True, range=(0.05,1.05))
        shuffled=True
    ax[0].plot(np.arange(len(cluster)),[0.8 for i in range(len(cluster))], visible=True, ls='--')
    ax[0].set_xlabel("cluster", fontsize=16)
    ax[0].set_ylabel("maximum fraction\nwith same %s"%label, fontsize=16)
    ax[0].set_ylim((0,1.1))
    ax[1].set_xlabel("maximum fraction\nwith same %s"%label, fontsize=16)
    ax[1].set_ylabel("pdf", fontsize=16)
    plt.show()
    fig.savefig("%s/%scluster_maximum_l%d_%s.pdf"%(directory,"shuffled" if shuffled else '',level,label))



def plot_maximum_size(clustersinfo, label, level, directory,clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(15,6))
    x = np.array(clustersinfo['sizes']).astype(int)
    y = np.array(clustersinfo['maximum'])[:,0].astype(float)
    plt.scatter(x,y, lw=10, label='clusters')
    plt.xlim(0,np.max(x)+np.max(x)/10)
    plt.plot(np.linspace(0.5,x.max()),1./np.linspace(0.5,x.max()), label='uniform')
    shuffled=False
    if clustersinfo_shuffle is not None:
        shuffled=True
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:,0].astype(float)
        plt.scatter(x_shuffle,y_shuffle, lw=10, label='clusters shuffled')
        plt.xlim(0,np.max(x_shuffle)+np.max(x_shuffle)/10)
    plt.xlabel("cluster size", fontsize=16)
    plt.ylabel("maximum fraction\nwith same %s"%label, fontsize=16)
    plt.ylim((0,1.1))
    plt.legend(loc='best', fontsize=18)
    plt.show()
    fig.savefig("%s/%sclusterhomosize_l%d_%s.pdf"%(directory, "shuffled" if shuffled else '', level,label))

def plot_maximum_label(clustersinfo,label,level, directory,clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(10,6))
    x = np.array(clustersinfo['nclasses']).astype(int)
    y = np.array(clustersinfo['maximum'])[:,0].astype(float)
    shuffled=False
    plt.scatter(x,y, lw=10, alpha=0.9, label='clusters')
    plt.plot(np.arange(1,np.max(x)+2),1./np.arange(1,np.max(x)+2),ls='--',c='cyan', label='uniform')
    plt.xlim(0.95,np.max(x)+0.5)
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:,0].astype(float)
        plt.scatter(x_shuffle,y_shuffle, lw=10, alpha=0.9, label='clusters shuffled')
        plt.plot(np.arange(1,np.max(x_shuffle)+2),1./np.arange(1,np.max(x_shuffle)+2),ls='--',c='cyan', label='')
        shuffled=True
        plt.xlim(0.95,np.max(x_shuffle)+0.5)
    plt.xlabel("# labels in cluster", fontsize=16)
    plt.ylabel("maximum fraction\nwith same %s"%label, fontsize=16)
    plt.ylim((0,1.1))
    plt.legend(loc='lower right', fontsize=18)
    plt.show()
    fig.savefig("%s/%scluster_homon_l%d_%s.pdf"%(directory, "shuffled" if shuffled else '', level,label))

def plot_labels_size(clustersinfo, label, level,directory, clustersinfo_shuffle=None):
    fig=plt.figure(figsize=(10,6))
    x = np.array(clustersinfo['sizes']).astype(float)
    y = np.array(clustersinfo['nclasses']).astype(int)
    plt.xlim(x.min()-10,x.max()+5)
    plt.ylim(y.min()-2,y.max()+5)
    shuffled=False
    plt.scatter(x,y, lw=10, alpha=0.9, label='clusters')
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(float)
        y_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        plt.scatter(x_shuffle,y_shuffle, lw=10, alpha=0.9, label='clusters shuffled')
        plt.xlim(x.min()-10,x_shuffle.max()+5)
        plt.ylim(y.min()-2,y_shuffle.max()+8)
        shuffled=True
    plt.xlabel("cluster size", fontsize=16)
    plt.ylabel("# labels in cluster", fontsize=16)
    plt.legend(loc='upper right', fontsize=18)
    plt.show()
    fig.savefig("%s/%scluster_shuffle_label_size_l%d_%s.pdf"%(directory, "shuffled" if shuffled else '', level,label))

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
    fig.savefig("%s/%sheatmap_cluster%s_l%d_%s.pdf"%(directory,"shuffled" if shuffled else '',"fraction_" if normalise else '', int(level),label))

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

def add_score_lines(ax, scores, labels, xl, h=False, c=False, alpha=0.2, **kwargs):
    '''
    add to ax lines in scores
    add homogeneity and completness if required by h and c
    '''
    colors = {
        'primary_site':'blue',
        'secondary_site':'red',
        'status':'red',
        'mixed':'green',
        'disease_type':'red',
        'shuffle': 'orange',
        'uniq':'purple',
        'hierarchical':'darkcyan',
        'lda':'violet',
        'RPPA Clusters':'red'
    }
    for label in labels:
        if h:
            ax.plot(xl, scores[label]['h'], ls='-.', c=colors[label], alpha=alpha, label='homogeneity - %s'%label)
        if c:
            ax.plot(xl, scores[label]['c'], ls=':', c=colors[label], alpha=alpha, label='completness - %s'%label)
        ax.plot(xl, scores[label]['V'], label='MI - %s'%label, ls='-', c=colors[label], **kwargs)
    customize_metric_plot(ax,xl)
        
def customize_metric_plot(ax, xl):
    ax.set_xlabel("number of clusters", fontsize=16)
    ax.set_ylabel("score", fontsize=16)
    ax.set_ylim((0,1.1))
    ax.set_xlim(np.min(xl),np.max(xl))
    ax.set_xscale('log')
    ax.legend(loc='best', fontsize=14)

def get_tissue_style(tissue):
    marker = 'o'
    c='k'
    ls='--'
    if 'gtex' in tissue:
        marker='o'
        ls='-'
    elif 'tcga' in tissue:
        marker='x'
        ls='--'
    else:
        marker='.'
        ls='-.'
    if 'reast' in tissue:
        c='darkcyan'
    elif 'olon' in tissue:
        c='b'
    elif 'hyroid' in tissue:
        c='y'
    elif 'terus' in tissue:
        c='pink'
    elif 'ladder' in tissue:
        c='gray'
    elif 'sophagus' in tissue:
        c='brown'
    elif 'ung' in tissue:
        c='magenta'
    elif 'tomach' in tissue:
        c='lime'
    elif 'kin' in tissue:
        c='wheat'
    elif 'ancreas' in tissue:
        c='forestgreen'
    elif 'Adrenal Gland' in tissue:
        c='aqua'
    elif 'Adipose Tissue' in tissue:
        c='brown'
    elif 'erve' in tissue:
        c='royalblue'
    elif 'lood' in tissue:
        c='red'
    elif 'idney' in tissue:
        c='mediumslateblue'
    elif 'eart' in tissue:
        c='darkred'
    elif 'rain' in tissue:
        c='darkgray'
    elif 'estis' in tissue:
        c='darkkhaki'
    else:
        c='k'
    return (marker,c,ls)