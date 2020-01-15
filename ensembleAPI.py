import mygene

mg = mygene.MyGeneInfo()

def geneinfo(genename):
    """Query ensemble for a ENSG gene

    Keyword arguments:
    genename -- ENSG
    """
    #gene = 123
    #genename = df['gene'][gene]
    try:
        q = mg.getgenes(genename[:15], 'name,symbol,refseq.rna,type_of_gene,bp')[0]
        print("%s [%s]"%(genename, q['symbol']))
        print("Descr: %s"%q['name'])
    except:
        pass

def genesinfo(genenames):
    """Query ensemble for a ENSG gene's array

    Keyword arguments:
    genename -- array like ENSG
    """
    #gene = 123
    #genename = df['gene'][gene]
    try:
        q_many = mg.getgenes(genenames, 'name,symbol,refseq.rna,type_of_gene,bp')
        for q in q_many:
            print("%s [%s]"%(q['query'], q['symbol']))
            print("Descr: %s"%q['name'])
    except:
        pass
