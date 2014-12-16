# -*- coding: utf-8 -*-

import os,sys,types
import Base, Genome
home = os.path.expanduser("~")
#datapath = os.path.join(home,'epitopedata')
datapath = os.path.join(request.folder,'static/data/')
genomespath = os.path.join(request.folder,'static/data/genomes')

def doTask():
    import time
    t = time.ctime()
    return dict(result=t)

def runPredictor(label,genome,names=None,methods='tepitope',
                alleles=['HLA-DRB1*0101'], **kwargs):
    """Run predictors"""

    if names != None:
        names = names.split(',')
    print label,genome,methods
    if type(methods) is not types.ListType:
        methods = [methods]
    for method in methods:
        P = Base.getPredictor(method)
        savepath = os.path.join(datapath, label)
        #if not os.path.exists(savepath):
        #    os.mkdir(savepath)
        gfile = os.path.join(genomespath,'%s.gb' %genome)
        df = Genome.genbank2Dataframe(gfile, cds=True)
        P.predictProteins(df,length=11,names=names,alleles=alleles,
                                label=label,save=True,path=savepath)
        #also pre-calculate binders
        b = Analysis.getAllBinders(savepath,method=method,n=n)
        b.to_csv(binderfile)
    rows = len(P.data)
    addPredictionstoDB(label)
    return dict(proteins=len(P.getNames()), rows=rows)

def addPredictionstoDB(label):
    """After a run is done add the prediction id to the db if not present"""

    return

def getOrthologs(seq):
    df = Analysis.getOrthologs(seq)
    return df

from gluon.scheduler import Scheduler
scheduler = Scheduler(db)

