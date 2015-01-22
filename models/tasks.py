# -*- coding: utf-8 -*-

import os,sys,types,glob
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
    #print row
    seq = row.translation
    prod = row['product']
    rec = SeqRecord(Seq(seq),id=tag,description=prod)
    fastafmt = rec.format("fasta")
    #print fastafmt
    feature = row.to_dict()
    ind = keys.get_loc(tag)
    previous = keys[ind-1]
    if ind<len(keys)-1:
        next = keys[ind+1]
    else:
        next=None
    return feature, fastafmt, previous, next

def getPredictions(label,genome,tag,q=0.96):
    """Get predictions from file system"""

    q=round(q,2)
    path = os.path.join(datapath, label)
    print path
    genomename = os.path.splitext(genome)[0]
    preds = OrderedDict()
    cutoffs = {}
    for m in methods:
        rpath = os.path.join(path, '%s/%s' %(genomename,m))
        filename = os.path.join(rpath, tag+'.mpk')
        print filename
        if not os.path.exists(filename):
            continue
        df = pd.read_msgpack(filename)
        pred = Base.getPredictor(name=m, data=df)
        #l=pred.getLength()
        cutoffs[m] = pred.allelecutoffs = Analysis.getCutoffs(rpath, m, q)
        if pred == None:
            continue
        preds[m] = pred
    return preds, cutoffs

def runPredictor(label,genome,newlabel='',names='',methods='tepitope',length=11,
                 mhc1alleles=[], drballeles=[], dpqalleles=[],
                 iedbmethod='IEDB_recommended',  **kwargs):
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
        P = Base.getPredictor(method)
        if method in ['iedbmhc1']:
            alleles = mhc1alleles
            P.iedbmethod = iedbmethod
            #P.path = iedbmhc1path
        else:
            alleles = drballeles + dpqalleles
        savepath = os.path.join(datapath, label, genome, method)
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        else:
            removeFiles(savepath)
        P.predictProteins(df,length=length,names=None,alleles=alleles,
                              label=label,save=True,path=savepath)
        #also pre-calculate binders for n=3
        b = Analysis.getAllBinders(savepath,method=method,n=3)
        binderfile = os.path.join(savepath, 'binders_3.csv')
        b.to_csv(binderfile)

    addPredictionstoDB(label,path=os.path.join(datapath, label))
    return dict(binders=len(b))

def removeFiles(path):
    filelist = glob.glob('*.mpk')
    for f in filelist:
        os.remove(f)
    return

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

def genomeAnalysis(label, gname, method, n=3):

    path = os.path.join(datapath, '%s/%s/%s' %(label,gname,method))

    if not os.path.exists(path):
       return dict(res=None)
    binderfile = os.path.join(path,'binders_%s.csv' %n)
    print binderfile
    if os.path.exists(binderfile):
        b = pd.read_csv(binderfile)
    else:
        b = Analysis.getAllBinders(path,method=method,n=n)
        b.to_csv(binderfile)

    P = Base.getPredictor(method)
    res = b.groupby('name').agg({P.scorekey:[np.mean,np.size,np.max]}).sort()
    res.columns = res.columns.get_level_values(1)

    #get genome and combine data
    gfile = getGenome(gname)
    g = Genome.genbank2Dataframe(gfile, cds=True)
    res = res.merge(g[['locus_tag','length','gene','product']],
                            left_index=True,right_on='locus_tag')
    #add protein urls
    res['locus_tag'] = res['locus_tag'].apply(
        lambda x: str(A(x, _href=URL(r=request,f='protein',args=[label,gname,x],extension=''))))
    res['perc'] = res.size/res.length*100
    res = res.sort('perc',ascending=False)
    print res[:10]
    cl = Analysis.findClusters(b, method)

    #get top binders in genome? plus most frequent
    #binders = b.groupby('peptide').agg({P.scorekey:np.size})
    #seabornsetup()
    fig=plt.figure()
    ax=fig.add_subplot(211)
    b.hist(P.scorekey,bins=30,alpha=0.8,ax=ax)
    ax.set_title('binder score distr')
    ax=fig.add_subplot(212)
    ax.set_title('coverage distr')
    res.hist('length',bins=30,alpha=0.8,ax=ax)
    plt.tight_layout()
    return b,res,cl,fig

from gluon.scheduler import Scheduler
scheduler = Scheduler(db)

