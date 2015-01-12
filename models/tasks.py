# -*- coding: utf-8 -*-

import os,sys,types
import ConfigParser
import Base, Genome, Analysis
home = os.path.expanduser("~")
datapath = os.path.join(request.folder,'static/results')
genomespath = os.path.join(request.folder,'static/data/genomes')

def doTask():
    import time
    t = time.ctime()
    return dict(result=t)

def runPredictor(label,genome,newlabel='',names='',methods='tepitope',
                 mhc1alleles=[], mhc2alleles=[],**kwargs):
    """Run predictors and save results"""

    applySettings()
    if names != ' ':
        names = names.split(',')
    if newlabel != '':
        label = newlabel
    query = db.genomes(db.genomes.name==genome)
    f,gfile = db.genomes.file.retrieve(query.file)
    if type(methods) is not types.ListType:
        methods = [methods]
    for method in methods:
        if method in ['iedbmhc2']:
            alleles = mhc1alleles
        else:
            alleles = mhc2alleles
        P = Base.getPredictor(method)
        savepath = os.path.join(datapath, label, method)
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        df = Genome.genbank2Dataframe(gfile, cds=True)
        P.predictProteins(df,length=11,names=None,alleles=alleles,
                                label=label,save=True,path=savepath)

        #also pre-calculate binders for n=3
        b = Analysis.getAllBinders(savepath,method=method,n=3)
        binderfile = os.path.join(savepath, 'binders_3.csv')
        b.to_csv(binderfile)

    addPredictionstoDB(label,path=os.path.join(datapath, label))
    return dict(binders=len(b))

def addPredictionstoDB(label,path):
    """Add the prediction id to the db if not present"""
    db.predictions.insert(identifier=label,path=path,user='')
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

def getConfig():
    """Get config file settings"""
    cpath = os.path.join(request.folder,'static')
    conffile = os.path.join(cpath, 'settings.conf')
    parser = ConfigParser.ConfigParser()
    parser.read(conffile)
    return parser, conffile

def applySettings():
    parser,conffile = getConfig()
    paths = dict(parser.items('base'))
    for i in paths:
        os.environ["PATH"] += os.pathsep+paths[i]
    return

from gluon.scheduler import Scheduler
scheduler = Scheduler(db)

