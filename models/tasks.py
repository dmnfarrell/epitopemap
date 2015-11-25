# -*- coding: utf-8 -*-

import os,sys,types,glob
import ConfigParser
from applications.epitopemap.modules.mhcpredict import base, sequtils, analysis
import pandas as pd
home = os.path.expanduser("~")
datapath = os.path.join(request.folder,'static/results')

def doTask():
    import time
    t = time.ctime()
    return dict(result=t)

def createMailer():
    """Create mailer for sending mails"""

    from gluon.tools import Mail
    #read in config
    parser,conffile = getConfig()
    settings = dict(parser.items('mail'))

    mail = Mail()
    mail.settings.server = settings['server']
    mail.settings.sender = settings['sender']
    mail.settings.login = None
    return mail

def sendMail(recipient, jobid):
    mail = createMailer()
    print mail
    mail.send(to=[recipient],
              subject='job %s finished' %jobid,
              message='hi there')

    return

def getGenome(name):
    """Get a genome file from the db"""

    record = db.genomes(db.genomes.name==name)
    filename = os.path.join(request.folder,'uploads',record.file)
    return filename

def getFasta(name):
    """Get a fasta file from db"""

    query = db.sequences(db.sequences.name==name)
    f,ffile = db.sequences.file.retrieve(query.file)
    desc = query.description
    return ffile, desc

def getTagbyGene(g,gene):
    """Get the locus_tag from the gene name"""

    fname = getGenome(g)
    df = sequtils.genbank2Dataframe(fname, cds=True)
    df = df.drop_duplicates('locus_tag')
    df = df.set_index('locus_tag')
    #if gene name provided we try to override the locus_tag
    t = df[df.gene==gene].index
    if len(t)>0:
        return t[0]
    else:
        return

def doSearch(genome, gene, desc):
    """Search genbank frame by gene or descr"""

    gfile = getGenome(genome)
    g = sequtils.genbank2Dataframe(gfile, cds=True)
    g = g.fillna('')
    if gene == '' and desc == '':
        df = g
    else:
        if gene == '': gene = '@#!x'
        if desc == '': desc = '@#!x'
        df = g[(g.gene.str.contains(gene, case=False)) |
                (g['product'].str.contains(desc, case=False))]
    df = df.drop(['type','pseudo','note','translation'],1)
    df = df.set_index('locus_tag')
    return df

def getFeature(g,tag):
    """Get gene feature from stored genbank file"""

    if g == 'other':
        return None,None,None,None
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    fname = getGenome(g)
    df = sequtils.genbank2Dataframe(fname, cds=True)
    df = df.drop_duplicates('locus_tag')
    df = df.set_index('locus_tag')
    keys = df.index
    if tag not in keys:
        return
    row = df.ix[tag]
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
    bcell = None
    for m in methods:
        rpath = os.path.join(path, '%s/%s' %(genomename,m))
        filename = os.path.join(rpath, tag+'.mpk')
        if not os.path.exists(filename):
            continue
        df = pd.read_msgpack(filename)
        pred = base.getPredictor(name=m, data=df)
        if m == 'bcell':
            bcell = pred
            continue
        cutoffs[m] = pred.allelecutoffs = analysis.getCutoffs(rpath, m, q)
        preds[m] = pred
    return preds, bcell, cutoffs

def runPredictors(label,genome='',newlabel='',names='',fasta='',methods='tepitope',length=11,
                 mhc1alleles=[], drballeles=[], dpqalleles=[], mhc2ref=False,
                 bcellmethod='Bepipred', iedbmethod='IEDB_recommended',
                 user=None, **kwargs):
    """Run predictors and save results"""

    limit = 40 #limit alleles
    applySettings()
    if names != '':
        names = names.split(',')
    else:
        names=None
    if newlabel != '':
        label = newlabel
    if genome != '':
        query = db.genomes(db.genomes.name==genome)
        f,gfile = db.genomes.file.retrieve(query.file)
        df = sequtils.genbank2Dataframe(gfile, cds=True)
    elif fasta != '':
        #if no genome given use the fasta sequences
        query = db.sequences(db.sequences.name==fasta)
        f,ffile = db.sequences.file.retrieve(query.file)
        df = sequtils.fasta2Dataframe(ffile)
        genome = 'other'
    else:
        return
    if type(methods) is not types.ListType:
        methods = [methods]
    length=int(length)
    for method in methods:
        P = base.getPredictor(method)
        if method in ['iedbmhc1']:
            alleles = mhc1alleles
            P.iedbmethod = iedbmethod
            #P.path = iedbmhc1path
        elif method == 'bcell':
            P.iedbmethod = bcellmethod
            alleles = None
        else:
            if type(drballeles) is types.StringType:
                drballeles=[drballeles]
            if type(dpqalleles) is types.StringType:
                dpqalleles=[dpqalleles]
            alleles = drballeles + dpqalleles
            if mhc2ref == 'on':
                alleles = base.getMHCIIRef()
        if len(alleles)>limit:
            alleles=alleles[:limit]
        savepath = os.path.join(datapath, label, genome, method)
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        else:
            removeFiles(savepath)
        P.predictProteins(df,length=length,names=names,alleles=alleles,
                              label=label,save=True,path=savepath)
        #also pre-calculate binders for n=3
        b = analysis.getAllBinders(savepath,method=method,n=3)
        if b is not None:
            binderfile = os.path.join(savepath, 'binders_3.csv')
            b.to_csv(binderfile)
    addPredictionstoDB(label,os.path.join(datapath, label),user)
    return dict(proteins=len(df))

def removeFiles(path):
    filelist = glob.glob('*.mpk')
    for f in filelist:
        os.remove(f)
    return

def addPredictionstoDB(label,path,user=''):
    """Add the prediction id to the db if not present"""

    qry = (db.predictions.identifier == label)
    rows = db(qry).select()
    if len(rows) == 0:
        db.predictions.insert(identifier=label,description='',user=user)
        db.commit()

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

def getBinders(path,method,n=3,cutoff=0.98):
    """Re-usable method to get binders from path"""

    binderfile = os.path.join(path,'binders_%s.csv' %n)
    print binderfile
    if os.path.exists(binderfile):
        b = pd.read_csv(binderfile)
    else:
        b = analysis.getAllBinders(path,method=method,n=n,cutoff=cutoff)
        b.to_csv(binderfile)
    return b

def genomeAnalysis(label, gname, method, n=3, cutoff=0.96):
    """Analyse over genome for cluster densities and other metrics"""

    path = os.path.join(datapath, '%s/%s/%s' %(label,gname,method))
    if not os.path.exists(path):
       return dict(res=None)
    b = getBinders(path,method,n,cutoff)
    P = base.getPredictor(method)
    if P.rankascending==1:
        func = np.min
    else:
        func = np.max
    res = b.groupby('name').agg({P.scorekey:[np.size,func]}).sort()
    res.columns = res.columns.get_level_values(1)
    res.rename(columns={'size': 'binders'}, inplace=True)

    #get genome and combine data
    gfile = getGenome(gname)
    g = sequtils.genbank2Dataframe(gfile, cds=True)
    res = res.merge(g[['locus_tag','length','gene','product','order']],
                            left_index=True,right_on='locus_tag')
    res['perc_binders'] = res['binders']/res.length*100
    res = res.sort('perc_binders',ascending=False)
    clusterfile = os.path.join(path,'clusters_%s.csv' %n)
    if os.path.exists(clusterfile):
        cl = pd.read_csv(clusterfile)
    else:
        cl = analysis.findClusters(b, method, dist=9)
        cl.set_index('name').to_csv(clusterfile)

    if cl is not None:
        gc = cl.groupby('name').agg({'density':np.max,'name':np.size})
        gc.rename(columns={'name': 'clusters'}, inplace=True)
        res = res.merge(gc,left_on='locus_tag',right_index=True,how='left')
    res = res.fillna('-')
    #get binder metrics
    top = b.groupby('peptide').agg({P.scorekey:func,'allele':np.max,
                    'name': lambda x: ','.join(list(x))}).reset_index()
    top = top.sort(P.scorekey,ascending=P.rankascending)[:1000]
    #add protein urls to results table
    res['locus_tag'] = res['locus_tag'].apply(
        lambda x: str(A(x, _href=URL(r=request,f='protein',args=[label,gname,x],extension=''))))
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
    preds, bcell, cutoffs = getPredictions(label,genome,tag,q=0.96)
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
        if genome != 'other':
            gfile = getGenome(genome)
            g = sequtils.genbank2Dataframe(gfile, cds=True)
            prot = g[g['locus_tag']==tag]
        else:
            ffile,desc = getFasta(genome)
            g = sequtils.fasta2Dataframe(ffile)
            prot = g[g['locus_tag']==tag]

        if len(prot)==0:
            return dict(res=None)
        seq = prot.translation.head(1).squeeze()
        print tag,seq
        alnrows = analysis.getOrthologs(seq,hitlist_size=400,equery=equery)
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

    fig=plt.figure(figsize=(5,5))
    ax=fig.add_subplot(111)
    stats=stats.T.drop_duplicates().T
    print stats
    stats.plot(ax=ax,grid=False,lw=1.5,colormap='rainbow')
    plt.legend(loc=2,prop={'size':9})
    ax.set_ylabel('percentage cons.')
    ax.set_title('conservation vs identity')
    plt.tight_layout()
    #fig = plotConserved(stats)
    return res, alnrows, summary, fig

def plotConserved(stats):
    from bokeh.models import HoverTool,ColumnDataSource
    from bokeh.plotting import Figure
    width=500
    height=500
    plot = Figure(title='',title_text_font_size="11pt",
            plot_width=width, plot_height=height,
           x_axis_label='identity',y_axis_label='perc_cons',
           tools="pan, wheel_zoom, resize, hover, reset, save",
           background_fill="#FAFAFA")
    x=stats.index
    for c in stats.columns:
        print x, stats[c]
        plot.line(x,stats[c], line_color='blue',line_width=1.5,legend=c)
    js,html = embedPlot(plot)
    return html

def correlation(label, genome, method1, method2, n=1):
    """Correlate binders from 2 methods"""

    n=int(n)
    path1 = os.path.join(datapath, '%s/%s/%s' %(label,genome,method1))
    path2 = os.path.join(datapath, '%s/%s/%s' %(label,genome,method2))
    print path1, path2
    if not os.path.exists(path1) or not os.path.exists(path2):
       return 1
    b1 = getBinders(path1,method1,n=n)
    b2 = getBinders(path2,method2,n=n)

    #m = pd.merge(b1, b2, on=['peptide','name','pos'])#, suffixes=['1','2']
    P1 = base.getPredictor(method1)
    P2 = base.getPredictor(method2)
    g1 = b1.groupby('name').agg({P1.scorekey:[np.mean,np.size,np.max]}).sort()
    g1.columns = g1.columns.get_level_values(1)
    g2 = b2.groupby('name').agg({P2.scorekey:[np.mean,np.size,np.max]}).sort()
    g2.columns = g2.columns.get_level_values(1)
    res = pd.merge(g1, g2, left_index=True,right_index=True)
    gfile = getGenome(genome)
    gn = sequtils.genbank2Dataframe(gfile, cds=True)
    res = res.merge(gn[['locus_tag','length']],
                            left_index=True,right_on='locus_tag')
    res['perc_x'] = res['size_x']/res.length*100
    res['perc_y'] = res['size_y']/res.length*100
    #print res[:20]
    return res

def findEpitope():
    """Find an epitope inside a set of seqs"""

    return

from gluon.scheduler import Scheduler
scheduler = Scheduler(db)
