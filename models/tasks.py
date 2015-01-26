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
    if tag not in keys:
        return
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
        #print filename
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

def runPredictors(label,genome,newlabel='',names='',methods='tepitope',length=11,
                 mhc1alleles=[], drballeles=[], dpqalleles=[],
                 iedbmethod='IEDB_recommended', **kwargs):
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

def genomeAnalysis(label, gname, method, n=3, cutoff=0.96):

    path = os.path.join(datapath, '%s/%s/%s' %(label,gname,method))

    if not os.path.exists(path):
       return dict(res=None)
    binderfile = os.path.join(path,'binders_%s.csv' %n)
    print binderfile
    if os.path.exists(binderfile):
        b = pd.read_csv(binderfile)
    else:
        b = Analysis.getAllBinders(path,method=method,n=n,cutoff=cutoff)
        b.to_csv(binderfile)

    P = Base.getPredictor(method)
    res = b.groupby('name').agg({P.scorekey:[np.mean,np.size,np.max]}).sort()
    res.columns = res.columns.get_level_values(1)

    #get genome and combine data
    gfile = getGenome(gname)
    g = Genome.genbank2Dataframe(gfile, cds=True)
    res = res.merge(g[['locus_tag','length','gene','product','order']],
                            left_index=True,right_on='locus_tag')
    res['perc'] = res['size']/res.length*100
    res = res.sort('perc',ascending=False)
    clusterfile = os.path.join(path,'clusters_%s.csv' %n)
    if os.path.exists(clusterfile):
        cl = pd.read_csv(clusterfile)
    else:
        cl = Analysis.findClusters(b, method, dist=9)
        cl.set_index('name').to_csv(clusterfile)

    if cl is not None:
        gc = cl.groupby('name').agg({'density':np.max})
        res = res.merge(gc,left_on='locus_tag',right_index=True,how='left')
    res = res.fillna('-')
    #get top binders in genome? plus most frequent
    top = b.groupby('peptide').agg({P.scorekey:np.max,'allele':np.max,
                    'name': lambda x: ','.join(list(x))}).reset_index()
    top = top.sort(P.scorekey,ascending=P.rankascending)[:1000]
    #add protein urls to results table
    res['locus_tag'] = res['locus_tag'].apply(
        lambda x: str(A(x, _href=URL(r=request,f='protein',args=[label,gname,x],extension=''))))
    #seabornsetup()
    #fig=plt.figure(figsize=(6,6))
    #ax=fig.add_subplot(111)
    #res.sort('order').plot(kind='bar',x='order',y='perc',ax=ax)
    #ax.set_title('binder score distr')
    #ax=fig.add_subplot(212)
    #ax.set_title('coverage distr')
    #res.hist('length',bins=30,alpha=0.8,ax=ax)
    #plt.tight_layout()
    fig=None
    return b,res,top,cl,fig

def conservationAnalysis(label, genome, method, tag, identity, n=3,
                         equery=None, **kwargs):
    """Conservation analysis"""

    n=int(n)
    identity=int(identity)
    res = {}
    blastpath = os.path.join(request.folder, 'static/data/blast')
    cachedfile = os.path.join(blastpath, '%s_%s.csv' %(genome,tag))

    #get predictions
    preds, cutoffs = getPredictions(label,genome,tag,q=0.97)
    if not preds.has_key(method):
        return 1
    pred = preds[method]
    length = pred.getLength()
    pb = pred.getPromiscuousBinders(n=n)
    if len(pb)==0:
        return 1
    if os.path.exists(cachedfile):
        print 'using %s' %cachedfile
        alnrows = pd.read_csv(cachedfile,index_col=0)
    else:
        gfile = getGenome(genome)
        g = Genome.genbank2Dataframe(gfile, cds=True)
        prot = g[g['locus_tag']==tag]
        if len(prot)==0:
            return dict(res=None)
        seq = prot.translation.head(1).squeeze()
        print tag,seq
        alnrows = Analysis.getOrthologs(seq,hitlist_size=400,equery=equery)
        if alnrows is None:
            alnrows=None
        #cache blast results for re-use
        alnrows.to_csv(cachedfile)

    alnrows.drop_duplicates(subset=['sequence'], inplace=True)

    def getConserved(pb,alnrows):
        #find conserved binders
        f=[]
        for i,a in alnrows.iterrows():
            seq = a.sequence
            found = [seq.find(j) for j in pb.peptide]
            f.append(found)
        s = pd.DataFrame(f,columns=pb.peptide,index=alnrows.accession)
        s = s.replace(-1,np.nan)
        res = s.count()
        return res

    def findConservedwithIdentity(alnrows,pb):
        vals=[]
        lowest = alnrows.perc_ident.min().round()
        for i in np.arange(lowest,100,10):
            x = alnrows[alnrows['perc_ident']>=i]
            s = getConserved(pb,x)/len(x)
            s.name=i
            vals.append(s)
        df = pd.DataFrame(vals)
        print df
        return df

    stats = findConservedwithIdentity(alnrows,pb)
    #limit to current identity level
    alnrows = alnrows[alnrows['perc_ident']>=identity]
    if len(alnrows) == 0:
        return 2
    found = getConserved(pb,alnrows)

    pb = pb.set_index('peptide')
    pb['conserved'] = found
    pb['perc_cons'] = pb.conserved/len(alnrows)
    pb=pb.sort('conserved',ascending=False).drop(['core','name'],1)
    summary = pb
    #print pb[:10]

    #seabornsetup()
    fig=plt.figure(figsize=(5,5))
    ax=fig.add_subplot(111)
    stats.plot(ax=ax,grid=False,lw=1.5)
    plt.legend(loc=2,prop={'size':9})
    ax.set_ylabel('percentage cons.')
    ax.set_title('conservation vs identity')
    plt.tight_layout()
    return res, alnrows, summary, fig


from gluon.scheduler import Scheduler
scheduler = Scheduler(db)

