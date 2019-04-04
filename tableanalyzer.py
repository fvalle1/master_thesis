import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def geneinfo(genename, df, nfiles):
    """Extimate mean and var for a ENSG gene

    Keyword arguments:
    genename -- ENSG
    genedata -- list or array-like object across tissues

    Return:
    genedict
    """
    #gene = 123
    #genename = df['gene'][gene]
    print("name: %s"%genename)
    genedata = np.array([fpkm for fpkm in df.loc[df['gene']==genename].loc[:,df.keys()[1:]].values.reshape(nfiles,1) if (fpkm>1e-1)&(fpkm<1e5)])
    try:
        genemean = np.nanmean(genedata)
        genevariance = np.nanvar(genedata)
        genedict = {
            'name':genename,
            'avg' : genemean,
            'var': genevariance,
            'data' : genedata,
        }
        try:
            q = mg.getgenes(genename[:15], 'name,symbol,refseq.rna,type_of_gene,bp')[0]
            print("Descr: %s"%q['name'])
            print("Symbol: %s"%q['symbol'])
            genedict['type']=q['type_of_gene']
        except:
            genedict['type']='unknown'
            pass
        print("FPKM mean: %10.2f"%genemean)
        print("FPKM var: %10.2f"%genevariance)
    except:
        return {}
    return genedict

def genedistr(genedict, bins = 50, ax = None, density=False, label='', save=True):
    """
    Plot distriution across tissues
    """
    maxfpkm = np.max(genedict['data'])
    width = float(maxfpkm) / bins
    _range = (0 - 0.5 * width, maxfpkm + 0.5 * width)
    if ax == None:
        fig = plt.figure(figsize=(15, 5))
        ax = fig.subplots()
    else:
        fig=ax.get_figure()
    n, bin_edges, _ = ax.hist(genedict['data'], lw=1.5, density=density, histtype='step', range=_range, bins=bins, label=label)
    ax.set_title(genedict['name'], fontsize=16)
    ax.set_xlabel('FPKM', fontsize=16)
    ax.set_ylabel('#', fontsize=16)
    ax.set_yscale('log')
    ax.set_xscale('log')
    if save:
        fig.savefig("plot/genes/%s_distr.pdf"%(genedict['name']))

def geneplot(genedict):
    """
    Plot FPKM across tissues
    """
    fig = plt.figure(figsize=(15, 5))
    plt.plot(genedict['data'], 'ob')
    plt.title(genedict['name'], fontsize=16)
    plt.xlabel("sample", fontsize=16)
    plt.ylabel("FPKM", fontsize=16)
    plt.yscale('log')
    plt.ylim(ymin=1e-4)
    plt.show()
    fig.savefig("plot/genes/%s_data.pdf"%(genedict['name']))

def genecoord(genedict, means, variances):
    """
    plot gene position in gobal plot
    """
    fig = plt.figure(figsize=(18,8))
    plt.scatter(means, variances)
    plt.scatter([np.average(genedict['data'])],[np.var(genedict['data'])], marker='x', c='r', s=90, label=genedict['name'])
    plt.xlabel("$<FPKM>$", fontsize=16)
    plt.ylabel("$\sigma^2_{FPKM}$", fontsize=16)
    plt.yscale('log')
    #plt.xlim(1e-3,200)
    plt.ylim(ymin=1e-2)
    plt.legend()
    plt.show()
    fig.savefig("plot/genes/%s_coord.png"%(genedict['name']))
