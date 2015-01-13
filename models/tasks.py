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

def getGenome(name):
    """Get a genome file from the db"""

    record = db.genomes(db.genomes.name==name)
    filename = os.path.join(request.folder,'uploads',record.file)
    return filename

def getFeature(g,tag):
    """Get gene feature from stored genbank file"""

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    fname = getGenome(g)
    df = Genome.genbank2Dataframe(fname, cds=True)
    df = df.drop_duplicates('locus_tag')
    df = df.set_index('locus_tag')
    keys = df.index
    row = df.ix[tag]
    print row
    seq = row.translation
    prod = row['product']
    rec = SeqRecord(Seq(seq),id=tag,description=prod)
    fastafmt = rec.format("fasta")
    print fastafmt
    feature = row.to_dict()
    ind = keys.get_loc(tag)
    previous = keys[ind-1]
    if ind<len(keys)-1:
        next = keys[ind+1]
    else:
        next=None
    return feature, fastafmt, previous, next

def getPredictions(label,genome,tag):
    """Get predictions from file system"""

    path = os.path.join(datapath, label)
    print path
    genomename = os.path.splitext(genome)[0]
    preds = OrderedDict()
    for m in methods:
        rpath = os.path.join(path, '%s/%s' %(genomename,m))
        filename = os.path.join(rpath, tag+'.mpk')
        print filename
        if not os.path.exists(filename):
            continue
        df = pd.read_msgpack(filename)
        pred = Base.getPredictor(name=m, data=df)
        l=pred.getLength()
        if m in Base.globalcutoffs:
            pred.allelecutoffs = Base.globalcutoffs[m][l]
        if pred == None:
            continue
        preds[m] = pred
    return preds

def runPredictor(label,genome,newlabel='',names='',methods='tepitope',length=11,
                 mhc1alleles=[], mhc2alleles=[],**kwargs):
    """Run predictors and save results"""

    applySettings()
    if names != ' ':
        names = names.split(',')
    if newlabel != '':
        label = newlabel
    query = db.genomes(db.genomes.name==genome)
    f,gfile = db.genomes.file.retrieve(query.file)
    df = Genome.genbank2Dataframe(gfile, cds=True)
    if type(methods) is not types.ListType:
        methods = [methods]
    length=int(length)
    for method in methods:
        if method in ['iedbmhc2']:
            alleles = mhc1alleles
        else:
            alleles = mhc2alleles
        P = Base.getPredictor(method)
        savepath = os.path.join(datapath, label, genome, method)
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        P.predictProteins(df,length=length,names=None,alleles=alleles,
                                label=label,save=True,path=savepath)

        #also pre-calculate binders for n=3
        b = Analysis.getAllBinders(savepath,method=method,n=3)
        binderfile = os.path.join(savepath, 'binders_3.csv')
        b.to_csv(binderfile)

    addPredictionstoDB(label,path=os.path.join(datapath, label))
    return dict(binders=len(b))

def addPredictionstoDB(label,path):
    """Add the prediction id to the db if not present"""
    db.predictions.insert(identifier=label,description='',user='')
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
    """Add binaries to path if needed"""
    parser,conffile = getConfig()
    paths = dict(parser.items('base'))
    for i in paths:
        os.environ["PATH"] += os.pathsep+paths[i]
    return

from gluon.scheduler import Scheduler
scheduler = Scheduler(db)

