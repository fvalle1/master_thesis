import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys
from mpl_finance import candlestick2_ohlc
import seaborn as sns
from sklearn import metrics


def plot_cluster_composition(fraction_sites, directory, level, normalise=False, label='primary_site', shuffled=False,
                             algorithm='topsbm'):
    sns.set(font_scale=0.8)
    df_clusters = pd.read_csv("%s/%s/%s_level_%d_clusters.csv" % (directory, algorithm, algorithm, level), header=[0])
    x = np.arange(1, 1 + len(df_clusters.columns))
    bottom = np.zeros(len(x))
    ymax = 0
    fig = plt.figure(figsize=(15, 8))
    for site, data in fraction_sites.items():
        if np.max(data) == 0:
            continue
        plt.bar(x, data, label=site, bottom=bottom)
        bottom = bottom + data
    plt.xlabel("cluster", fontsize=20)
    if normalise:
        plt.ylabel("fraction of nodes", fontsize=20)
    else:
        plt.ylabel("number of nodes", fontsize=20)
    plt.title("%s%s distribution across clusters" % ("Shuffled " if shuffled else '', label), fontsize=20)
    plt.legend(ncol=3)
    plt.xticks(x)
    plt.show()
    fig.savefig("%s/%s/%s%sclustercomposition_l%d_%s.pdf" % (
        directory, algorithm, "shuffled" if shuffled else '', "fraction_" if normalise else '', int(level), label))


def get_cluster_given_l(l, directory, algorithm='topsbm'):
    df_clusters = pd.read_csv("%s/%s/%s_level_%d_clusters.csv" % (directory, algorithm, algorithm, l), header=[0],
                              index_col=None)
    cluster = {}
    for i, c in enumerate(df_clusters.columns):
        cluster[i] = df_clusters[c].dropna().values
    return cluster


def get_topic_given_l(l, directory, algorithm='topsbm'):
    df_topics = pd.read_csv("%s/%s/%s_level_%d_topics.csv" % (directory, algorithm, algorithm, l), header=[0])
    topic = {}
    for i, c in enumerate(df_topics.columns):
        topic[i] = df_topics[c].dropna().values
    return topic


def get_fraction_sites(cluster, df_files, label='primary_site', normalise=False):
    fraction_sites = {}
    c_fraction_site = {}
    for site in df_files[label].dropna().unique():
        fraction_sites[site] = []
        c_fraction_site[site] = 0

    for i, c in enumerate(cluster):
        for sample in cluster[i]:
            try:
                for fullsample in df_files.index.values:
                    if sample in fullsample:
                        foundsample = df_files.loc[fullsample, :]
                c_fraction_site[foundsample[label]] += 1
            except:
                print("error in %s" % sample)
        for site in fraction_sites:
            if normalise:
                norm = float(len(cluster[i]))
            else:
                norm = 1
            fraction_sites[site].append(c_fraction_site[site] / norm)
            c_fraction_site[site] = 0
    df = pd.DataFrame(data=fraction_sites)
    ##put first columns that have high values in average
    avgs = df.apply(lambda x: np.average(x.to_numpy()[x.to_numpy().nonzero()[0]]), axis=0)
    df = df.transpose()
    df.insert(0, 'avg', avgs)
    df = df.sort_values(by=['avg'], axis=0, ascending=False).drop('avg', axis=1).transpose()
    df = df.sort_values(by=[tissue for tissue in df.columns], axis=0, ascending=False)
    return df.to_dict(orient='list')


def get_clustersinfo(cluster, fraction_sites):
    clustersinfo = {
        "maximum": [],
        "homogeneity": [],
        "sizes": [],
        "nclasses": []
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
                nclass += 1
                # using fraction_items normalised
                if cdata <= 1:
                    homo -= cdata * np.log(cdata)
            size += cdata
        if cumulative > 0:
            clustersinfo['maximum'].append([float(maximum) / cumulative, site_maximum])
        else:
            clustersinfo['maximum'].append([0, site_maximum])
        clustersinfo['sizes'].append(size)
        clustersinfo['nclasses'].append(nclass)
        clustersinfo['homogeneity'].append(1 - homo)
    return clustersinfo


def plot_maximum(clustersinfo, cluster, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(15, 6))
    ax = fig.subplots(1, 2)
    bins = 10
    real = np.array(clustersinfo['maximum'])[:, 0].astype(float)
    ax[0].plot(np.sort(real), marker='o', ms=25, ls='')
    ax[1].hist(np.sort(real), histtype='step', bins=bins, lw=4, density=True, range=(0.05, 1.05))
    shuffled = False
    if clustersinfo_shuffle is not None:
        shuffled = np.array(clustersinfo_shuffle['maximum'])[:, 0].astype(float)
        ax[0].plot(np.sort(shuffled), marker='o', ls='', ms=25)
        ax[1].hist(np.sort(shuffled), histtype='step', bins=bins, lw=4, density=True, range=(0.05, 1.05))
        shuffled = True
    ax[0].plot(np.arange(len(cluster)), [0.8 for i in range(len(cluster))], visible=True, ls='--')
    ax[0].set_xlabel("cluster", fontsize=20)
    ax[0].set_ylabel("maximum fraction\nwith same %s" % label, fontsize=20)
    ax[0].set_ylim((0, 1.1))
    ax[1].set_xlabel("maximum fraction\nwith same %s" % label, fontsize=20)
    ax[1].set_ylabel("pdf", fontsize=20)
    plt.show()
    fig.savefig(
        "%s/%s/%scluster_maximum_l%d_%s.pdf" % (directory, algorithm, "shuffled" if shuffled else '', level, label))


def plot_maximum_size(clustersinfo, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(15, 6))
    x = np.array(clustersinfo['sizes']).astype(int)
    y = np.array(clustersinfo['maximum'])[:, 0].astype(float)
    plt.scatter(x, y, lw=10, label='clusters')
    plt.xlim(0, np.max(x) + np.max(x) / 10)
    plt.plot(np.linspace(0.5, x.max()), 1. / np.linspace(0.5, x.max()), label='uniform')
    shuffled = False
    if clustersinfo_shuffle is not None:
        shuffled = True
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:, 0].astype(float)
        plt.scatter(x_shuffle, y_shuffle, lw=10, label='clusters shuffled')
        plt.xlim(0, np.max(x_shuffle) + np.max(x_shuffle) / 10)
    plt.xlabel("cluster size", fontsize=20)
    plt.ylabel("maximum fraction\nwith same %s" % label, fontsize=20)
    plt.ylim((0, 1.1))
    plt.legend(loc='best', fontsize=20)
    plt.show()
    fig.savefig(
        "%s/%s/%sclusterhomosize_l%d_%s.pdf" % (directory, algorithm, "shuffled" if shuffled else '', level, label))


def plot_maximum_label(clustersinfo, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(10, 6))
    x = np.array(clustersinfo['nclasses']).astype(int)
    y = np.array(clustersinfo['maximum'])[:, 0].astype(float)
    shuffled = False
    plt.scatter(x, y, lw=10, alpha=0.9, label='clusters')
    plt.plot(np.arange(1, np.max(x) + 2), 1. / np.arange(1, np.max(x) + 2), ls='--', c='cyan', label='uniform')
    plt.xlim(0.95, np.max(x) + 0.5)
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:, 0].astype(float)
        plt.scatter(x_shuffle, y_shuffle, lw=10, alpha=0.9, label='clusters shuffled')
        plt.plot(np.arange(1, np.max(x_shuffle) + 2), 1. / np.arange(1, np.max(x_shuffle) + 2), ls='--', c='cyan',
                 label='')
        shuffled = True
        plt.xlim(0.95, np.max(x_shuffle) + 0.5)
    plt.xlabel("number of labels", fontsize=20)
    plt.ylabel("maximum fraction\nwith same %s" % label, fontsize=20)
    plt.ylim((0, 1.1))
    plt.legend(loc='lower right', fontsize=20)
    plt.show()
    fig.savefig(
        "%s/%s/%scluster_homon_l%d_%s.pdf" % (directory, algorithm, "shuffled" if shuffled else '', level, label))


def plot_labels_size(clustersinfo, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(10, 6))
    x = np.array(clustersinfo['sizes']).astype(float)
    y = np.array(clustersinfo['nclasses']).astype(int)
    plt.xlim(x.min() - 10, x.max() + 5)
    plt.ylim(y.min() - 2, y.max() + 5)
    shuffled = False
    plt.scatter(x, y, lw=10, alpha=0.9, label='clusters')
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(float)
        y_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        plt.scatter(x_shuffle, y_shuffle, lw=10, alpha=0.9, label='clusters shuffled')
        plt.xlim(x.min() - 10, x_shuffle.max() + 5)
        plt.ylim(y.min() - 2, y_shuffle.max() + 8)
        shuffled = True
    plt.xlabel("cluster size", fontsize=20)
    plt.ylabel("number of labels", fontsize=20)
    plt.legend(loc='upper right', fontsize=20)
    plt.show()
    fig.savefig(
        "%s/%s/%scluster_shuffle_label_size_l%d_%s.pdf" % (
            directory, algorithm, "shuffled" if shuffled else '', level, label))


def make_heatmap(fraction_sites, directory, label, level, shuffled=False, normalise=False, algorithm='topsbm'):
    sns.set(font_scale=2)
    found_classes = []
    for site, data in fraction_sites.items():
        if np.max(data) == 0:
            continue
        found_classes.append(site)
    for arr in fraction_sites.values():
        x = len(arr)
        break
    x = np.arange(1, 1 + x)
    fig = plt.figure(figsize=(30, 10))
    fig.subplots(1)
    sns.heatmap(pd.DataFrame(data=fraction_sites).loc[:, found_classes].transpose(), vmin=0, cmap="RdYlBu_r",
                xticklabels=x)
    fig.savefig("%s/%s/%sheatmap_cluster%s_l%d_%s.pdf" % (
        directory, algorithm, "shuffled" if shuffled else '', "fraction_" if normalise else '', int(level), label))


def get_file(sample, df_file):
    for fullsample in df_file.index.values:
        if sample in fullsample:
            return df_file.loc[fullsample, :]


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
                print(sys.exc_info()[0])
                print("error searching %s in %s" % (label,sample))
    _, true_labels = np.unique(true_labels, return_inverse=True)
    return true_labels, predicted_labels


def add_score_lines(ax, scores, labels=None, xl=[], h=False, c=False, alpha=0.2, **kwargs):
    '''
    add to ax lines in scores
    add homogeneity and completness if required by h and c
    '''
    colors = {
        'primary_site': 'blue',
        'secondary_site': 'red',
        'status': 'red',
        'hSBM': 'green',
        'disease_type': 'red',
        'shuffle': 'orange',
        'disease_tissue': 'purple',
        'hierarchical': 'darkcyan',
        'lda': 'violet',
        'RPPA Clusters': 'red'
    }
    for label in labels:
        if h:
            ax.plot(xl, scores[label]['h'], ls='-.', c=colors[label], alpha=alpha, label='homogeneity - %s' % label)
        if c:
            ax.plot(xl, scores[label]['c'], ls=':', c=colors[label], alpha=alpha, label='completness - %s' % label)
        if len(scores[label]['V']) == len(xl):
            ax.plot(xl, scores[label]['V'], label='MI - %s' % label, ls='-', c=colors[label], **kwargs)
    customize_metric_plot(ax, xl)


def customize_metric_plot(ax, xl):
    ax.set_xlabel("number of clusters", fontsize=20)
    ax.set_ylabel("V-measure score", fontsize=20)
    ax.set_ylim((0, 1.1))
    ax.set_xlim(np.min(xl), np.max(xl))
    ax.set_xscale('log')
    ax.legend(loc='best', fontsize=20)


def plot_topic_size(directory, l, algorithm='topsbm'):
    df_topics = pd.read_csv("%s/%s/%s_level_%d_topics.csv" % (directory, algorithm, algorithm, l))
    sizes = []
    for t in df_topics.columns:
        sizes.append(len(df_topics.loc[:, t].dropna()))
    bins = np.linspace(0.5, np.max(sizes) + 0.5, int((np.max(sizes) + 1) / (np.max(sizes) / 5)))
    bin_counts, bin_edges, _ = plt.hist(sizes, histtype='step', lw=2, bins=bins)
    fig = plt.figure()
    ax = fig.subplots()
    ax.set_title("[%d topics, level: %d]" % (len(df_topics.columns), l))
    x = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax.plot(x[np.nonzero(bin_counts)], bin_counts[np.nonzero(bin_counts)])
    ax.plot(x, 1e4 / np.power(x, 5))
    ax.set_xlabel("topic size\n(number of genes)", fontsize=20)
    ax.set_ylabel("number of topic", fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()
    fig.savefig("%s/%s/topic_size_level%d.png" % (directory, algorithm, l))


def get_candles(directory, level, df_mv, ax, algorithm='topsbm'):
    df_topics = pd.read_csv("%s/%s/%s_level_%d_topics.csv" % (directory, algorithm, algorithm, level))
    candles = {
        'open': [],
        'high': [],
        'low': [],
        'close': [],
        'size': []
    }
    for topic in df_topics.columns:
        subarr = df_mv.loc[df_topics[topic].dropna(), :]['occurrence'].values
        avg = np.average(subarr)
        std = np.std(subarr)
        q = np.quantile(subarr, [0.25, 0.75])
        candles['high'].append(np.min([1, avg + std]))
        candles['open'].append(np.min([q[1], 1]))
        candles['close'].append(np.max([q[0], 0]))
        candles['low'].append(np.max([0, avg - std]))
        candles['size'].append(len(subarr))
    ax.set_title("[level: %d]" % level)
    ax.set_ylabel('$O_i$', fontsize=20)
    ax.set_xlim(-1, len(df_topics.columns))
    ax.set_xticks([i + 1 for i in range(-1, len(df_topics.columns))])
    ax.set_xticklabels(
        ["Topic %d" % (i + 2) if ((i + 2) % 5 == 0 or i == -1) else '' for i in range(-1, len(df_topics.columns))],
        rotation=60)
    return candles


def get_tissue_style(tissue):
    marker = 'o'
    c = 'k'
    ls = '--'
    if 'gtex' in tissue:
        marker = 'o'
        ls = '-'
    elif 'tcga' in tissue:
        marker = 'x'
        ls = '--'
    else:
        marker = '.'
        ls = '-.'
    if 'reast' in tissue:
        c = 'darkcyan'
    elif 'olon' in tissue:
        c = 'b'
    elif 'hyroid' in tissue:
        c = 'y'
    elif 'terus' in tissue:
        c = 'pink'
    elif 'ladder' in tissue:
        c = 'gray'
    elif 'sophagus' in tissue:
        c = 'brown'
    elif 'ung' in tissue:
        c = 'magenta'
    elif 'tomach' in tissue:
        c = 'lime'
    elif 'kin' in tissue:
        c = 'wheat'
    elif 'ancreas' in tissue:
        c = 'forestgreen'
    elif 'Adrenal Gland' in tissue:
        c = 'aqua'
    elif 'Adipose Tissue' in tissue:
        c = 'brown'
    elif 'erve' in tissue:
        c = 'royalblue'
    elif 'lood' in tissue:
        c = 'red'
    elif 'idney' in tissue:
        c = 'mediumslateblue'
    elif 'eart' in tissue:
        c = 'darkred'
    elif 'rain' in tissue:
        c = 'darkgray'
    elif 'estis' in tissue:
        c = 'darkkhaki'
    else:
        c = 'k'
    return (marker, c, ls)


def topic_distr_sample(doc, df, ax=None):
    if ax == None:
        fig = plt.figure()
        ax = fig.subplots()
    ax.set_title("Topic distribution: %s" % doc)
    labels = [l if df[df['doc'] == doc].loc[:, l].values[0] >= 0.05 else '' for l in df.columns[2:]]
    patches, texts, autotexts = ax.pie(df[df['doc'] == doc].values[0][2:], labels=labels,
                                       autopct=lambda p: '%.1f%s' % (p, '%') if p >= 5 else '',
                                       textprops={'fontsize': 20, 'color': 'white', 'wrap': True})
    for t in texts:
        t.set_fontsize(18)
        t.set_wrap(True)
        t.set_color('black')
    plt.show()


def topic_distr_isample(idoc, df, ax=None):
    topic_distr_sample(df[df['i_doc'] == idoc]['doc'].values[0], ax)

def add_tumor_location(df_files):
    df_files.insert(2, 'disease_tissue', '')
    for sample in df_files.index.values:
        row = df_files.loc[sample, :]
        df_files.at[sample, 'disease_tissue'] = '%s[%s]' % (row['primary_site'], row['disease_type'])

def get_scores(directory, labels, L=3, verbose=False):
    df_files = pd.read_csv("%s/files.dat" % directory, index_col=[0], header=[0])
    if df_files.columns.isin(['disease_type']).any():
        add_tumor_location(df_files)
    scores = {}
    for label in labels:
        scores[label] = {
            'h': [],
            'c': [],
            'V': []
        }
        for l in np.arange(L + 1):
            true_labels, predicted_labels = define_labels(get_cluster_given_l(l, directory), df_files, label=label)
            scores[label]['h'].append(metrics.cluster.homogeneity_score(true_labels, predicted_labels))
            scores[label]['c'].append(metrics.cluster.completeness_score(true_labels, predicted_labels))
            scores[label]['V'].append(metrics.cluster.v_measure_score(true_labels, predicted_labels))
            try:
                if verbose:
                    print(l)
            except:
                pass
    if len(labels) >= 2:
        h = np.array(scores[labels[0]]['h'])
        c = np.array(scores[labels[1]]['c'])
        scores['mixed'] = {
            'h': h,
            'c': c,
            'V': 2 * h * c / (h + c)
        }
    scores['shuffle'] = {
        'h': [],
        'c': [],
        'V': []
    }
    try:
        for l in np.arange(0, L + 1):
            if verbose:
                print(l)
            _, predicted_labels = define_labels(get_cluster_given_l(l, directory), df_files, label='primary_site')
            true_labels, _ = define_labels(get_cluster_given_l(l, directory),
                                           pd.read_csv("%s/files.dat.shuf" % directory, index_col=[0]),
                                           label='primary_site')
            scores['shuffle']['h'].append(metrics.cluster.homogeneity_score(true_labels, predicted_labels))
            scores['shuffle']['c'].append(metrics.cluster.completeness_score(true_labels, predicted_labels))
            scores['shuffle']['V'].append(metrics.cluster.v_measure_score(true_labels, predicted_labels))
    except:
        print("shuffled files not found")
    return scores


def getclustersizesarray(directory, L=3, algorithm='topsbm'):
    xl = []
    try:
        xl = [len(get_cluster_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(0, L, L + 1)]
    except:
        try:
            xl = [len(get_cluster_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(1, L, L)]
        except:
            raise ("error saving clustersizes")
    return xl


def gettopicsizesarray(directory, l=3, algorithm='topsbm'):
    xl = []
    try:
        xl = [len(get_topic_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(0, l, l + 1)]
    except:
        try:
            xl = [len(get_topic_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(1, l, l)]
        except:
            raise ("error saving clustersizes")
    return xl


def clusteranalysis(directory, labels, l=3, algorithm='topsbm'):
    df_clusters = pd.read_csv("%s/%s/%s_level_%d_clusters.csv" % (directory, algorithm, algorithm, l), header=[0])
    if df_clusters.isnull():
        print("files not found")
    df_files = pd.read_csv("%s/files.dat" % directory, index_col=[0], header=[0])
    for normalise in [True, False]:
        for label in labels:
            for level in np.arange(l + 1)[::-1]:
                if level == 0:
                    continue
                print(normalise, label, level)
                try:
                    cluster = get_cluster_given_l(level, directory, algorithm=algorithm)
                    fraction_sites = get_fraction_sites(cluster, df_files=df_files, label=label, normalise=normalise)

                    # fsdf = pd.DataFrame(data=fraction_sites)
                    # fsdf = fsdf.drop('Other', axis=1)
                    # fsdf = fsdf.divide(fsdf.sum(axis=1), axis=0).fillna(0)
                    # fraction_sites = fsdf.sort_values(by=fsdf.columns.to_list(), ascending=True).to_dict(orient='list')

                    clustersinfo = get_clustersinfo(cluster, fraction_sites)
                    plot_cluster_composition(fraction_sites, directory, level, label=label, normalise=normalise,
                                             algorithm=algorithm)
                    make_heatmap(fraction_sites, directory, label, level, normalise=normalise, algorithm=algorithm)

                    if not normalise:
                        plot_maximum(clustersinfo, cluster, label, level, directory, algorithm=algorithm)
                        plot_maximum_size(clustersinfo, label, level, directory, algorithm=algorithm)
                        plot_maximum_label(clustersinfo, label, level, directory, algorithm=algorithm)
                except:
                    print(sys.exc_info()[0])
                try:
                    fraction_sites_shuffle = get_fraction_sites(cluster,
                                                                pd.read_csv("%s/files.dat.shuf" % directory,
                                                                            index_col=[0]),
                                                                label=label, normalise=normalise)
                    clustersinfo_shuffle = get_clustersinfo(cluster, fraction_sites_shuffle)
                    plot_cluster_composition(fraction_sites_shuffle, directory, level, normalise=normalise, label=label,
                                             shuffled=True)
                    if not normalise:
                        plot_maximum(clustersinfo, cluster, label, level, directory, clustersinfo_shuffle,
                                     algorithm=algorithm)
                        plot_maximum_size(clustersinfo, label, level, directory, clustersinfo_shuffle, algorithm=algorithm)
                        plot_maximum_label(clustersinfo, label, level, directory, clustersinfo_shuffle, algorithm=algorithm)
                        plot_labels_size(clustersinfo, label, level, directory, clustersinfo_shuffle, algorithm=algorithm)
                except:
                    print("must shuffle files")

    ##define scores
    scores = get_scores(directory, labels)
    try:
        xl = getclustersizesarray(directory, l)
        with open("%s/clustersizes.txt" % directory, 'w') as f:
            for x in xl:
                f.write("%d\n" % x)
    except:
        print("cannot save clustersizes.txt")

    try:
        xl = gettopicsizesarray(directory, l)
        with open("%s/topicsizes.txt" % directory, 'w') as f:
            for x in xl:
                f.write("%d\n" % x)
    except:
        print("cannot save topicsizes.txt")

    # save files for R analisys
    for l in np.arange(l + 1):
        pd.DataFrame(data=define_labels(get_cluster_given_l(l, directory), df_files, label=labels[0])[1],
                     columns=['l%d' % l]).to_csv("%s/%s/%s_level_%d_labels.csv" % (directory, algorithm, algorithm, l),
                                                 header=True,
                                                 index=False)
