# -*- coding: utf-8 -*-

import os,sys,types
import Base, Genome, Analysis
home = os.path.expanduser("~")
#datapath = os.path.join(home,'epitopedata')
datapath = os.path.join(request.folder,'static/data/results')
genomespath = os.path.join(request.folder,'static/data/genomes')

def doTask():
    import time
    t = time.ctime()
    return dict(result=t)

def runPredictor(label,genome,newlabel=None,names=None,methods='tepitope',
                 mhc1alleles=[], mhc2alleles=[],**kwargs):
    """Run predictors and save results"""

    if names != None:
        names = names.split(',')
    if newlabel != None:
        label = newlabel
    query = db.genomes(db.genomes.name==genome)
    f,gfile = db.genomes.file.retrieve(query.file)
    #print gfile
    if type(methods) is not types.ListType:
        methods = [methods]
    for method in methods:
        if method in ['iedbmhc2']:
            alleles = mhc1alleles
        else:
            alleles = mhc2alleles
        P = Base.getPredictor(method)
        savepath = os.path.join(datapath, label, method)
        df = Genome.genbank2Dataframe(gfile, cds=True)
        P.predictProteins(df,length=11,names=None,alleles=alleles,
                                label=label,save=True,path=savepath)

        #also pre-calculate binders for n=3
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

def getSeqDepot(seq):
    """Fetch seqdepot annotation for sequence"""
    import SeqDepot
    sd = SeqDepot.new()
    aseqid = sd.aseqIdFromSequence(seq)
    try:
        result = sd.findOne(aseqid)
    except Exception, e:
        print e
        result=None
    return result

from gluon.scheduler import Scheduler
scheduler = Scheduler(db)

