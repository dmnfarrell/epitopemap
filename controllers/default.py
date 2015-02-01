# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## This is a samples controller
## - index is the default action of any application
## - user is required for authentication and authorization
## - download is for downloading files uploaded in the db (does streaming)
## - call exposes all registered services (none by default)
#########################################################################

from gluon.contrib.user_agent_parser import mobilize
import os,sys,types
import string,operator
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib as mpl
import bokeh
#from bokeh.plotting import *
home = os.path.expanduser("~")
genomespath = os.path.join(request.folder,'static/data/genomes')
#datapath = os.path.join(home,'epitopedata')
datapath = os.path.join(request.folder,'static/results')
import Genome, Base, Tepitope, Analysis
methods = ['tepitope','netmhciipan','iedbmhc1','bcell']#,'threading'] #'iedbmhc2'
iedbmethods = ['IEDB_recommended','consensus','ann','smm','arb','netmhcpan']
bcellmethods = ['Chou-Fasman', 'Emini', 'Karplus-Schulz',
                'Kolaskar-Tongaonkar', 'Parker', 'Bepipred']
colors = {'tepitope':'green','netmhciipan':'orange',
           'iedbmhc1':'blue','iedbmhc2':'pink','threading':'purple'}
colormaps={'tepitope':'Greens','netmhciipan':'Reds','iedbmhc2':'Oranges',
               'threading':'Purples','iedbmhc1':'Blues'}

def index():
    """
    example action using the internationalization operator T and flash
    rendered by views/default/index.html or views/generic.html
    """
    if request.user_agent().is_mobile:
        response.view.replace('.html','.mobile.html')
    form = quicksearch()
    return dict(message=T('Menu'),searchform=form)

def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    """
    auth.settings.registration_requires_approval = True
    return dict(form=auth())

def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request,db)

def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()

@auth.requires_signature()
def data():
    return dict(form=crud())

def mpld3Plot(fig, objects=None):
    """mpld3 html plot from figure"""

    import mpld3
    html = mpld3.fig_to_html(fig)
    htmllabels = []
    if objects!=None and len(objects)>0:
        bars,labels = zip(*objects)
        tooltip = MyToolTip(bars, labels)
        plugins.connect(fig, tooltip)
    return html

def embedPlot(plot):
    """Embed plot method for new versions of bokeh"""

    from bokeh.resources import Resources
    from bokeh.embed import autoload_static
    fp = os.path.join(request.folder,'static/temp/')
    fp = os.path.join(fp, plot._id+".js")
    res = Resources("relative")
    res.js_files = ["../static/js/bokeh.min.js"]
    res.css_files = ["../static/css/bokeh.min.css"]
    jspath = os.path.join('../static/temp/', plot._id+".js")
    js,tag = autoload_static(plot, res, jspath)
    with open(fp, "w") as f:
        f.write(js)
    print
    return js,tag

def plotAnnotations(plot,annotation):
    h=1.8
    y=.4+h/2.0
    if 'tmhmm' in annotation:
        vals = annotation['tmhmm']
        x=[i[0]+(i[1]-i[0])/2.0 for i in vals]
        w=[i[1]-i[0] for i in vals]
        print x,w,y
        #source = ColumnDataSource(data=dict(x=x,y=y))
        plot.rect(x,y, width=w, height=h,color='blue',line_color='blue',alpha=0.6)
    if 'pfam27' in annotation:
        vals = annotation['pfam27']
        print vals
        text = [i[0] for i in vals]
        x=[i[1]+(i[2]-i[1])/2.0 for i in vals]
        w=[i[2]-i[1] for i in vals]
        print x,w,y
        plot.rect(x,y, width=w, height=h,color='white',line_color='black',alpha=0.6)
        plot.text(x,y-1, text=text, angle=0)
    return

def plotBCell(plot,pred,height):
    """Line plot of b cell predictions - no allele stuff"""

    x = pred.data.Position
    #print pred.data[:20]
    #source = ColumnDataSource(data=dict(x=x,y=y))
    y=pred.data.Score
    h=height
    y = y+abs(min(y))
    y = y*(h/max(y))+3
    plot.line(x, y, line_color="red", line_width=2, alpha=0.6,legend='bcell')
    return

def plotLines(preds,tag,width=820,height=None,seqdepot=None,bcell=None):
    """Plot epitopes as lines"""
    from bokeh.models import Range1d,HoverTool,ColumnDataSource
    from bokeh.plotting import Figure

    yrange = Range1d(start=1, end=10)
    if height==None:
        height = 130+50*len(preds)
    plot = Figure(title=tag,title_text_font_size="11pt",plot_width=width, plot_height=height,
           y_range=yrange,
           tools="xpan, xwheel_zoom, resize, hover, reset, save",
           background_fill="#FAFAFA")
    for m in preds:
        pred = preds[m]
        sc = pred.scorekey
        y=pred.data[sc]

    return

def plotTracks(preds,tag,n=3,title=None,width=820,height=None,seqdepot=None,bcell=None):
    """Plot epitopes as parallel tracks"""

    from bokeh.models import Range1d,HoverTool,FactorRange,Grid,GridPlot,ColumnDataSource
    from bokeh.plotting import Figure

    alls=1
    if title == None:
        title=tag
    for m in preds:
        alls += len(preds[m].data.groupby('allele'))
    if height==None:
        height = 130+10*alls
    yrange = Range1d(start=0, end=alls+3)
    plot = Figure(title=title,title_text_font_size="11pt",plot_width=width, plot_height=height,
           y_range=yrange, #y_range=ylabels,
           y_axis_label='allele',
           tools="xpan, xwheel_zoom, resize, hover, reset, save",
           background_fill="#FAFAFA")

    h=3
    if bcell != None:
        plotBCell(plot, bcell, alls)
    if seqdepot != None:
        plotAnnotations(plot,seqdepot)

    #lists for hover data
    #we plot all rects at once
    x=[];y=[];allele=[];widths=[];clrs=[];peptide=[]
    predictor=[];position=[];score=[];leg=[]

    for m in preds:
        pred = preds[m]
        cmap = mpl.cm.get_cmap(colormaps[m])
        df = pred.data
        sckey = pred.scorekey
        pb = pred.getPromiscuousBinders(data=df,n=n)
        #print pb
        l = pred.getLength()
        grps = df.groupby('allele')
        alleles = grps.groups.keys()
        if len(pb)==0:
            continue
        c=colors[m]
        leg.append(m)

        for a,g in grps:
            b = pred.getBinders(data=g)
            b = b[b.pos.isin(pb.pos)] #only promiscuous
            b.sort('pos',inplace=True)
            scores = b[sckey].values
            score.extend(scores)
            pos = b['pos'].values
            position.extend(pos)
            x.extend(pos+(l/2.0)) #offset as coords are rect centers
            widths.extend([l for i in scores])
            clrs.extend([c for i in scores])
            y.extend([h+0.5 for i in scores])
            alls = [a for i in scores]
            allele.extend(alls)
            peptide.extend(list(b.peptide.values))
            predictor.extend([m for i in scores])
            h+=1

    source = ColumnDataSource(data=dict(x=x,y=y,allele=allele,peptide=peptide,
                                    predictor=predictor,position=position,score=score))
    plot.rect(x,y, width=widths, height=0.8,
         x_range=Range1d(start=1, end=len(g)+l),
         color=clrs,line_color='gray',alpha=0.7,source=source,
         min_border_left=2,min_border_right=2)

    hover = plot.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("allele", "@allele"),
        ("position", "@position"),
        ("peptide", "@peptide"),
        ("score", "@score"),
        ("predictor", "@predictor"),
    ])
    #seqlen = pred.data.pos.max()+l
    #plot.set(x_range=Range1d(start=1, end=seqlen))

    plot.xaxis.major_label_text_font_size = "8pt"
    plot.xaxis.major_label_text_font_style = "bold"
    plot.ygrid.grid_line_color = None
    plot.yaxis.major_label_text_font_size = '0pt'
    plot.xaxis.major_label_orientation = np.pi/4

    js,html = embedPlot(plot)
    #print html
    return html

def plotEmpty(width=850):
    """Plot an empty plot"""

    from bokeh.models import Range1d
    plot = figure(title='',plot_width=width, plot_height=10,
           y_range=Range1d(start=1, end=100),
           tools="xpan, xwheel_zoom, resize, hover, reset",
           background_fill="white")
    x=range(100); y=2
    rect(x,y, width=1, height=0.8,color='white')
    js,html = embedPlot(plot)
    print plot
    return html

def plots():
    """Use as component to plot predictions for given request"""

    print 'plot request'
    print request.vars
    label = request.vars.label
    #if we have no data
    if label == 'dummy':
        figure = plotEmpty()
        return dict(figure=figure)
    g = request.vars.genome
    tag = request.vars.tag
    gene = request.vars.gene
    title=None
    if gene != None:
        t = getTagbyGene(g,gene) #override tag with gene name if provided
        if t != None:
            tag = t
            title = tag+' / '+gene
    if request.vars.width == None:
        width = 820
    else:
        width = int(request.vars.width)
    if request.vars.height != None:
        height = int(request.vars.height)
    else:
        height = None
    if request.vars.n == None:
        n=3
    else:
        n = int(request.vars.n)
    if request.vars.perccutoff != None:
        perccutoff=float(request.vars.perccutoff)
    else:
        perccutoff=0.97
    preds,bcell,cutoffs = getPredictions(label,g,tag,perccutoff)
    if len(preds)==0 or preds==None:
        return dict(error=True)
    sd=None
    if request.vars.annotation == 'on':
        feat, fastafmt, previous, next = getFeature(g,tag)
        seq = feat['translation']
        sd = getSeqDepot(seq)['t']

    if request.vars.kind == None:
        figure = plotTracks(preds,tag,n=n,title=title,
                    width=width,height=height,seqdepot=sd,bcell=bcell)
    else:
        figure = plotLines(preds,tag,width=width,height=height,seqdepot=sd,bcell=bcell)

    return dict(figure=figure,preds=preds,error=False)

def scoredistplots(preds):
    """Score distribution plots"""

    plots=[]
    for p in preds:
        pred=preds[p]
        key=pred.scorekey
        data = pred.data[key]
        hist, edges = np.histogram(data, density=True, bins=30)
        p = figure(title=p,plot_height=250,tools='')
        p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
               fill_color="#036564", line_color="#033649")
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        plots.append(p)
    plot = GridPlot(children=[plots],title='test')
    js,html = embedPlot(plot)
    return html

def results():
    """Component to show predictions for all peptides for each predictor """

    label = request.vars.label
    g = request.vars.genome
    tag = request.vars.tag
    preds,bcell,cutoffs = getPredictions(label,g,tag)
    summary = summaryhtml(preds)
    data = {}
    for p in preds:
        data[p] = preds[p].reshape()
    data = dict(data)
    return dict(data=data)

def binders():
    """Component for top binder tables"""

    label = request.vars.label
    g = request.vars.genome
    tag = request.vars.tag
    n = int(request.vars.n)
    preds,bcell,cutoffs = getPredictions(label,g,tag)
    summary = summaryhtml(preds)
    b = Base.getBinders(preds,n=n)
    kys = b.keys()
    if 'tepitope' in kys and 'netmhciipan' in kys:
        shared = pd.merge(b['tepitope'],b['netmhciipan'],
                    on=['peptide','name','pos','core'],
                    copy=False).sort('pos')
    else:
        shared=''
    return dict(b=b,summary=summary,shared=shared,n=n)

def showSequence(seq,preds):
    """Get html display of binders on sequences"""

    colors = {'tepitope':'#70E2AA','netmhciipan':'#FF8181',
              'iedbmhc1':'#9DCEFF','iedbmhc2':'orange','threading':'#BCA9F5'}
    l=9 #need to get this from predictors
    seqs=[]
    tabledata=[]
    #idx =  ''.join([seq[i] if i%10!=0 else '|' for i in range(len(seq))])
    tabledata.append((TR(TH('allele'),TH('sequence'))))
    for p in preds:
        b = preds[p].getBinders()
        clr = colors[p]
        #pb = preds[p].getPromsicuousBinders(n=n)
        #b = b[b.pos.isin(pb.pos)]
        grps = b.groupby('allele')
        for a,g in grps:
            pos=[]
            for i in g.pos: pos.extend(np.arange(i,i+l))
            seqhtml=[]
            for i in range(len(seq)):
                if i in pos:
                    seqhtml.append(SPAN(seq[i],_style="background-color:%s" %clr))
                else:
                    seqhtml.append(SPAN(seq[i],_style="color: gray"))
            tabledata.append((TR(TH(a),TD(*seqhtml))))
    table = TABLE(*tabledata,_class="seqtable")
    return table

def sequence():
    """Component to highlight epitopes on sequence"""

    label = request.vars.label
    g = request.vars.genome
    tag = request.vars.tag
    n = int(request.vars.n)
    feat, fastafmt, previous, next = getFeature(g,tag)
    seq = feat.qualifiers['translation'][0]
    preds,bcell,c = getPredictions(label,g,tag)
    table = showSequence(seq,preds)
    return dict(table=table)

def feature():
    """Component showing gene annotation"""

    g = request.vars.genome
    tag = request.vars.tag
    items = getFeature(g,tag)
    if items != None:
        feat, fastafmt, previous, next = items
    return dict(fastafmt=fastafmt,feat=feat,
                previous=previous,next=next)
    return dict()

def iedb():
    """remote iedb tools predcitions"""

    g = request.vars.genome
    tag = request.vars.tag
    feature, fastafmt, previous, next = getFeature(g,tag)
    seq = feature.qualifiers['translation'][0]
    df = Base.getIEDBRequest(seq)
    result = XML(df.to_html(classes='mytable'))
    return dict(result=result)

def seqdepot(result):
    """Sedepot data table format"""

    #print result
    kys = result['t'].keys()
    tables = {}
    for k in kys:
        fieldnames = [TH(i) for i in sd.toolFields(k)]
        rows = [TR(i) for i in result['t'][k]]
        rows.insert(0,TR(*fieldnames))
        tables[k] = TABLE(*rows,_class="tinytable")
    fp = os.path.join(request.folder,'static/temp/')
    filename = os.path.join(fp,tag+'.png')
    sd.saveImage(aseqid, filename, {'format':'png'})
    imgurl = IMG(_src=URL(r=request,c='static',f='temp/%s' %os.path.basename(filename)))
    links = [LI(A(k,_href="#%s" %k)) for k in tables]
    tablinks = UL(*links,_class="small-tab-links")
    divs=[DIV(tables[k],_id=k,_class="tab") for k in tables]
    content = DIV(*divs,_class="tab-content")
    tabbedcontent = DIV(tablinks, content,_class="tabs")
    return dict(result=result,seq=seq,imgurl=imgurl,tables=tables,
                tabbedcontent=tabbedcontent)

def protein():
    """Display protein info from a fixed URL"""

    label = request.args[0]
    g = request.args[1]
    tag = request.args[2]
    n = 3

    items = getFeature(g,tag)
    if items != None:
        feature, fastafmt, previous, next = items
    else:
        raise HTTP(404, "No such feature %s available in genome %s" %(tag,g))
        return

    result = dict(label=label,tag=tag,genome=g,n=n,
                   previous=previous,next=next)
    return result

@auth.requires_login()
def genomes():
    """Display available genomes and allow upload"""

    formats = ['genbank'] #,'embl','fasta']
    uploadform = FORM(
                   TABLE(TR(TD(LABEL('Identifier:',_for='name')),
                        TD(INPUT(_name='name',_type='string',_required=True))),
                     TR(TD(LABEL('Format:',_for='format')),
                        TD(SELECT(formats,_name='format',_type='string',_required=True))),
                     TR(TD(LABEL('file to upload')),TD(INPUT(_name='gfile',_type='file'))),
                     TR(TD(),TD(INPUT(_name='submit',_type='submit',_value='Submit'))),
                     _class="smalltable"), _id="myform")

    if uploadform.accepts(request.vars,formname='upload_form'):
		fname = request.vars.gfile.filename
		uploadform.vars.filename = fname
		id = db.genomes.insert(name=uploadform.vars.name,
		                    file=uploadform.vars.gfile,
		                    filename=uploadform.vars.filename,
		                    format=uploadform.vars.format)

    db.genomes.id.readable=False
    query=((db.genomes.id>0))
    default_sort_order=[db.genomes.id]
    links=[lambda row: A('browse',_href=URL('genome', args=row.name))]
    grid = SQLFORM.grid(query=query,  orderby=default_sort_order,
                create=False, deletable=True, maxtextlength=350, paginate=35,
                details=True, csv=False, ondelete=myondelete,
                editable=auth.has_membership('editor_group'),links=links)

    return dict(grid=grid,form=uploadform)

def genome():
    """Summary page for genome"""

    g = request.args[0]
    if len(request.args) == 1:
        gfile = getGenome(g)
        data = Genome.genbank2Dataframe(gfile)
        summary = Genome.genbankSummary(data)
        data=data[data.type=='CDS']
        data=data.set_index('locus_tag')
        return dict(genome=g,data=data,summary=summary)
    else:
        return dict()

@auth.requires_login()
def predictions():
    """Manage available predictions"""

    '''
    if addform.validate():
        path = addform.vars.path
        addform.vars.identifier = os.path.basename(path)
        session.flash = 'form accepted'
        db.predictions.insert(**dict(addform.vars))
    if addform.process().accepted:
        session.flash = 'form accepted'
    '''
    db.predictions.id.readable=False
    query=((db.predictions.id>0))
    default_sort_order=[db.predictions.id]
    grid = SQLFORM.grid(query=query, orderby=default_sort_order,
                create=True, deletable=True, maxtextlength=350,
                paginate=20,details=True, csv=False, ondelete=myondelete,
                editable=auth.has_membership('editor_group'))

    return dict(grid=grid)

def myondelete(table, id):
    form = FORM.confirm('Are you sure?')
    print form
    if form.accepted:
        response.flash = "I don't like your submission"
        print table, id
        #db(db.predictions.id==id).delete()
    return form

def summaryhtml(predictors):
    """Summary table of predictions"""

    rows=[]
    rows.append(TR(TH('name'),TH('cutoff'),TH('binders')))
    for p in predictors:
        pred=predictors[p]
        b = pred.getPromiscuousBinders(n=2)
        rows.append(TR(pred.name, pred.cutoff, len(b)))
    return TABLE(*rows,_class='tinytable')

def download():

    import StringIO
    label = request.args[0]
    g = request.args[1]
    t = request.args[2]
    preds,bcell,c = getPredictions(label,g,t)
    data = [preds[p].data for p in preds]
    df = pd.concat(data)
    output = StringIO.StringIO()
    df.to_csv(output,float_format='%.2f')
    csvdata = output.getvalue()
    return dict(csvdata=csvdata)

def clusterResults():
    """Cluster results"""

    results = {}
    files = ['topclusters_MTB-H37Rv.csv','topsinglebinders.csv']
    for f in files:
        f = os.path.join(datapath, f)
        r = pd.read_csv(f, index_col=0)
        r.reset_index(inplace=True,drop=True)
        r.sort('name',inplace=True)
        results[f] = r
    return dict(results=results)

def search():
    """Advanced search form"""

    form = SQLFORM.factory(
              Field("genome", requires=IS_IN_DB(db, 'genomes.id', '%(name)s',
                    zero=None,multiple=False,orderby=~db.genomes.name)),
                  Field("locus_tag", default='Rv0011c'),
                  Field("gene", default=''),
                  Field("product", default='',comment='any words from product field'),
                  Field("start", type='integer', default='',comment='genomic coords'),
                  Field("end", type='integer', default=''),
                  formstyle='table3cols',
                  submit_button="Search")
    results={}
    msg = T(" ")
    # testing if the form was accepted
    if form.process().accepted:
        query = db.proteins.id>0
        # gathering form submitted values
        if form.vars.locus_tag :
            query &= db.proteins.locus_tag == form.vars.locus_tag
        if form.vars.gene:
            query &= db.proteins.gene == form.vars.gene
        if form.vars.type:
            query &= db.proteins.type == form.vars.type
        if form.vars.product:
            query &= db.proteins.product == form.vars.product
        if form.vars.start:
            query &= db.proteins.start >= form.vars.start
        if form.vars.end:
            query &= db.proteins.end <= form.vars.end

        count = db(query).count()
        results = db(query).select(orderby=~db.proteins.locus_tag)
        msg = T("%s records found" %count )

    return dict(form=form,msg=msg,results=results)

def quicksearch():
    """Non DB search just using paths"""

    form = SQLFORM.factory(
              Field('label',requires=IS_IN_DB(db, 'predictions.identifier',zero=None,
                    multiple=False),default=1,label='id'),
              Field('genome',requires=IS_IN_DB(db, 'genomes.name', zero=None,
                    multiple=False),default=1,label='genome'),
              Field('tag', 'string', label='locus tag',default='',length=10),
              Field('gene', 'string', label='gene',default='',length=10),
              hidden=dict(width=550,height=250,n=2),
              formstyle="table3cols",_id='myform')
    form.element('input[name=tag]')['_style'] = 'width:210px;'
    form.element('input[name=gene]')['_style'] = 'width:210px;'
    return form

def selectionForm():

    form = SQLFORM.factory(
              Field('label',requires=IS_IN_DB(db, 'predictions.identifier',zero=None,
                    multiple=False),default=1,label='id'),
              Field('genome',requires=IS_IN_DB(db, 'genomes.name', zero=None,
                    multiple=False),default=1,label='genome'),
              Field('tag', 'string', label='locus tag',default=''),
              Field('n', 'string', label='min alleles',default=3),
              Field('globalcutoff', 'boolean', label='global cutoff',default=True),
              Field('perccutoff', 'string', label='perc. cutoff',default=.98),
              Field('annotation', 'boolean', label='annotation',default=False),
              submit_button="Update",
              formstyle='table3cols',_id='myform',_class='myform')
    form.element('input[name=n]')['_style'] = 'width:50px;'
    form.element('input[name=perccutoff]')['_style'] = 'width:50px;'
    #form.element('input[name=scorecutoff]')['_style'] = 'width:50px;'
    form.element('input[name=tag]')['_style'] = 'width:130px;'
    return form


def quickview():
    """Quickview"""

    defaultid = 'results_bovine'
    form = selectionForm()
    return dict(label=defaultid,form=form)

def show():
    """Quickview all results in one - faster"""

    label = request.vars.label
    g = request.vars.genome
    tag = request.vars.tag
    n = int(request.vars.n)
    cutoff = float(request.vars.perccutoff)
    print request.vars
    if request.vars.perccutoff == None:
        cutoff = 0.95
    else:
        cutoff = float(request.vars.perccutoff)
    if request.vars.width == None:
        width = 820
    else:
        width = int(request.vars.width)
    annot = request.vars.annotation

    if label == 'dummy':
        figure = plotEmpty()
    preds,bcell,cutoffs = getPredictions(label,g,tag,cutoff)

    if len(preds)==0:
        redirect(URL('error'))

    feat, fastafmt, previous, next = getFeature(g,tag)
    seq = feat['translation']
    sd=None
    if request.vars.annotation == 'on':
        sd = getSeqDepot(seq)['t']
    figure = plotTracks(preds,tag,n=n,width=width,seqdepot=sd,bcell=bcell)
    #distplots = scoredistplots(preds)
    summary = summaryhtml(preds)
    #get all results into tables
    data = {}
    for p in preds:
        data[p] = preds[p].reshape()
    data = dict(data)
    #top binders
    b = Base.getBinders(preds,n=n)
    kys = b.keys()
    if 'tepitope' in kys and 'netmhciipan' in kys:
        shared = pd.merge(b['tepitope'],b['netmhciipan'],
                    on=['peptide','name','pos','core'],
                    copy=False).sort('pos')
    else:
        shared=''
    seqtable = showSequence(seq,preds)
    #info
    path = os.path.join(datapath, label)
    found = [(m,preds[m].getLength()) for m in preds]
    info = TABLE(*found,_class='tinytable')

    return dict(figure=figure,feat=feat,fastafmt=fastafmt,data=data,#distplots=distplots
                b=b,summary=summary,shared=shared,n=n,seqtable=seqtable,cutoffs=cutoffs,
                genome=g,tag=tag,label=label,info=info,path=path)

def error():
    return dict()

def analysis():
    """Genome wide analysis of epitope predictions"""
    defaultid = 'results_test'
    predids = [p.identifier for p in db().select(db.predictions.ALL)]
    opts1 = [OPTION(i,value=i) for i in predids]
    genomes = [p.name for p in db().select(db.genomes.ALL)]
    opts2 = [OPTION(i,value=i) for i in genomes]
    form = FORM(TABLE(
            TR(TD(LABEL('id:',_for='genome')),
            TD(SELECT(*opts1,_name='label',
                    value=defaultid, _style="width:150px;"))),
            TR(TD(LABEL('genome:',_for='genome')),
            TD(SELECT(*opts2,_name='genome',value='',_style="width:150px;"))),
            TR(TD(LABEL('method:',_for='method')),
            TD(SELECT(*methods,_name='method',value='tepitope',_style="width:150px;"))),
            TR(TD(LABEL('min alleles:',_for='n')),
            TD(INPUT(_name='n',_type='text',value=3,_style="width:50px;"))),
            TR(TD(LABEL('perc cutoff:',_for='perccutoff')),
            TD(INPUT(_name='perccutoff',_type='text',value='0.96',_style="width:50px;"))),
            TR(TD(),TD(INPUT(_name='submit',_type='submit',_value='Analyse'))),
            _class="smalltable"), _id="myform")
    return dict(form=form)

@auth.requires_login()
def analysegenome():
    """Analyse genome predictions"""

    pd.set_option('max_colwidth', 800)
    gname = request.vars.genome
    label = request.vars.label
    method = request.vars.method
    n = int(request.vars.n)
    cutoff = float(request.vars.perccutoff)

    b,res,top,cl,fig = genomeAnalysis(label, gname, method, n, cutoff)
    #plothtml = mpld3Plot(fig)
    plothtml=''
    summary = '%s proteins with %s binders in >%s alleles' %(len(res),len(b),n)
    return dict(res=res,top=top,cl=cl,summary=summary,plothtml=plothtml)

def compare():
    """Correlate predictions from 2 methods"""

    form = SQLFORM.factory(
              Field('label',requires=IS_IN_DB(db, 'predictions.identifier',zero=None,
                    multiple=False),default=1,label='id'),
              Field('genome',requires=IS_IN_DB(db, 'genomes.name', zero=None,
                    multiple=False),default=1,label='genome'),
              Field('method1',requires=IS_IN_SET(methods,multiple=False,zero=None),label='method 1'),
              Field('method2',requires=IS_IN_SET(methods,multiple=False,zero=None),label='method 2'),
              Field('n', 'string', label='min alleles',default=3),
              hidden=dict(perccutoff=.98),
              formstyle="table3cols",_id='myform',_class='myform')
    form.element('input[name=n]')['_style'] = 'width:50px;'
    return dict(form=form)

def correlationanalysis():

    fig=''
    msg=None
    if request.vars.method1 == request.vars.method2:
        return dict(res=None,msg='2 methods are the same!')
    print request.vars
    res = correlation(**request.vars)
    if res is None:
        msg = 'no such predictions'
    fig = plotCorrelation(res)
    return dict(fig=fig,res=res,msg=msg)

def plotCorrelation(res):
    from bokeh.models import HoverTool,ColumnDataSource
    from bokeh.plotting import Figure
    width=600
    height=600
    plot = Figure(title='',title_text_font_size="11pt",
            plot_width=width, plot_height=height,
           x_axis_label='method1',y_axis_label='method2',
           tools="pan, wheel_zoom, resize, hover, reset, save",
           background_fill="#FAFAFA")

    x=res['perc_x']
    y=res['perc_y']
    source = ColumnDataSource(data=dict(x=x,y=y, protein=res.locus_tag))
    plot.circle(x,y, color='blue', line_color='gray',fill_alpha=0.5, size=10, source=source)
    hover = plot.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("binders1", "@x"),
        ("binders2", "@y"),
        ("protein", "@protein"),
    ])
    js,html = embedPlot(plot)
    return html

def conservationAnalysisForm(defaultid='test'):

    defaultg = 'MTB-H37Rv'
    predids = [p.identifier for p in db().select(db.predictions.ALL)]
    opts1 = [OPTION(i,value=i) for i in predids]
    genomes = [p.name for p in db().select(db.genomes.ALL)]
    opts2 = [OPTION(i,value=i) for i in genomes]
    form = FORM(TABLE(
            TR(TD(LABEL('id:',_for='genome')),
            TD(SELECT(*opts1,_name='label',
                    value=defaultid, _style="width:150px;"))),
            TR(TD(LABEL('genome:',_for='genome')),
            TD(SELECT(*opts2,_name='genome',value=defaultg,_style="width:150px;"))),
            TR(TD(LABEL('locus tag:',_for='tag')),
            TD(INPUT(_name='tag',_type='text',value="Rv0001",_style="width:150px;"))),
            TR(TD(LABEL('method:',_for='method')),
            TD(SELECT(*methods,_name='method',value='tepitope',_style="width:150px;"))),
            TR(TD(LABEL('min alleles:',_for='n')),
            TD(INPUT(_name='n',_type='text',value=3,_style="width:50px;"))),
            TR(TD(LABEL('min identity:',_for='identity')),
            TD(INPUT(_name='identity',value=70,_style="width:50px;"))),
            TR(TD(),TD('BLAST options')),
            TR(TD(LABEL('entrez query:',_for='entrezquery')),
            TD(TEXTAREA(_name='entrezquery',value='',_style="height:100px;width:150px;"))),
            TR(TD(),TD(INPUT(_name='submit',_type='submit',_value='Submit'))),
            _class="smalltable"), _id="myform", hidden=dict(width=850))
    return form

def conservation():
    """Analysis of epitope conservation"""

    form = conservationAnalysisForm()
    '''if form.process().accepted:
        session.flash = 'form accepted'
        pvars = {'seq':seq,'hitlist_size':400,'equery':equery}
        task = scheduler.queue_task('doTask', #pvars=request.vars,
                                    immediate=True, timeout=300)
        print task.id
        status = scheduler.task_status(task.id, output=True)
        result = status.result
        print status'''
    return dict(form=form)

@auth.requires_login()
def conservationanalysis():
    """Analysis of epitope conservation"""

    pd.set_option('max_colwidth', 3000)
    label = request.vars.label
    gname = request.vars.genome
    method = request.vars.method
    n=int(request.vars.n)
    tag = request.vars.tag
    identity = int(request.vars.identity)
    equery = request.vars.entrezquery
    retval = conservationAnalysis(**request.vars)
    msg=''
    if retval == 1:
        msg =  'No predictions found for %s with method %s with n=%s.' %(tag,method,n)
        return dict(res=None,msg=msg)
    elif retval == 2:
        msg =  'No BLAST results at >%s%% sequence identity.' %identity
        return dict(res=None,msg=msg)
    else:
        res, alnrows, summary, fig = retval

    alnrows = Analysis.getAlignedBlastResults(alnrows)
    alnrows = Analysis.setBlastLink(alnrows)
    plothtml = mpld3Plot(fig)
    url =  A('direct link to these results', _href=URL('default','conservationanalysis.load',
                vars={'label':label,'genome':gname,'tag':tag,'method':method,'n':n,
                 'identity':identity,'equery':equery},extension=''))
    return dict(res=res,alnrows=alnrows,summary=summary,plothtml=plothtml,
                msg=msg,permlink=url)

def submissionForm():
    """Form for job submission"""

    applySettings() #so that paths to predictors work
    defaultg = 'MTB-H37Rv'
    predids = [p.identifier for p in db().select(db.predictions.ALL)]
    opts1 = [OPTION(i,value=i) for i in predids]
    genomes = [p.name for p in db().select(db.genomes.ALL)]
    opts2 = [OPTION(i,value=i) for i in genomes]
    p1 = Base.getPredictor('iedbmhc1')
    mhc1alleles = p1.getMHCIList()
    p2 = Base.getPredictor('netmhciipan')
    mhc2alleles = p2.getAlleleList()
    drballeles = Base.getDRBList(mhc2alleles)
    dqpalleles = Base.getDQPList(mhc2alleles)
    tepitopealleles = Tepitope.getAlleles()
    #get all possible alleles for both MHCII methods
    drballeles = sorted(list(set(drballeles+tepitopealleles)))
    lengths = [9,11,13,15]
    form = FORM(DIV(
            TABLE(
            TR(TD(LABEL('current labels:',_for='genome')),
            TD(SELECT(*opts1,_name='label',
                    value='', _style="width:200px;"))),
            TR(TD(LABEL('new label:',_for='genome')),
            TD(INPUT(_name='newlabel',_type='text',value="",_style="width:200px;"))),
            TR(TD(LABEL('genome:',_for='genome')),
            TD(SELECT(*opts2,_name='genome',value=defaultg,_style="width:200px;"))),
            TR(TD(LABEL('locus tags:',_for='names')),
            TD(INPUT(_name='names',_type='text',value="",_style="width:200px;"))),
            TR(TD(LABEL('methods:',_for='methods')),
            TD(SELECT(*methods,_name='methods',value='tepitope',_size=4,_style="width:200px;",
                _multiple=True))),
            TR(TD(LABEL('mhc1 method:',_for='iedbmethod')),
            TD(SELECT(*iedbmethods,_name='iedbmethod',value='IEDB_recommended',_size=1,
                _style="width:200px;"))),
            TR(TD(LABEL('bcell method:',_for='bcellmethod')),
            TD(SELECT(*bcellmethods,_name='bcellmethod',value='Bepipred',_size=1,
                _style="width:200px;"))),
            TR(TD(LABEL('length:',_for='length')),
            TD(SELECT(*lengths,_name='length',value=11,_size=1,_style="width:70px;"))),
            TR(TD(),TD(INPUT(_name='submit',_type='submit',_value='Submit Job'))),
            _class="smalltable"),_style='float: left'),
            DIV(TABLE(
            TR(TD(LABEL('MHC-I alleles:',_for='alleles')),
            TD(SELECT(*mhc1alleles,_name='mhc1alleles',value='HLA-A*01:01-10',_size=7,_style="width:200px;",
                _multiple=True))),
            TR(TD(LABEL('MHC-II DRB:',_for='alleles')),
            TD(SELECT(*drballeles,_name='drballeles',value='HLA-DRB1*0101',_size=7,_style="width:200px;",
                _multiple=True))),
            TR(TD(LABEL('MHC-II DQ/P:',_for='alleles')),
            TD(SELECT(*dqpalleles,_name='dqpalleles',value='',_size=6,_style="width:200px;",
                _multiple=True))),
            _class="smalltable"),_style='float: left'),

            _id="myform")#, _class='myform')

    return form

@auth.requires_login()
def submit():
    """Process job for submission and queue job"""

    form = submissionForm()
    if form.process().accepted:
        session.flash = 'form accepted'
        task = scheduler.queue_task('runPredictors', pvars=request.vars,
                                    immediate=True, timeout=86400)
        redirect(URL('jobsubmitted', vars={'id':task.id}))
    elif form.errors:
        response.flash = 'form has errors'
    return dict(form=form)

@auth.requires_login()
def jobsubmitted():
    """Get details of a submitted job"""

    taskid = int(request.vars['id'])
    status = scheduler.task_status(taskid, output=True)
    return dict(taskid=taskid,status=status)

def test():
    l='results_bovine'
    g='MTB-H37Rv'
    tag='Rv0011c'
    preds,bcell,c = getPredictions(l,g,tag)
    #html = plotTracks(preds,tag)
    form = FORM(TABLE(
            TD(INPUT(_name='submit',_type='submit',_value='Update'))), _id="myform")
    return dict(g=g,tag=tag,l=l,form=form)

def bokehtest():
    """Bokeh test"""

    from bokeh.models import Range1d, HoverTool, GridPlot
    from bokeh.plotting import Figure
    N = 100
    x = np.random.random(size=N) * 100
    y = np.random.random(size=N) * 100
    radii = np.random.random(size=N) * 3
    colors = ["#%02x%02x%02x" % (r, g, 150) for r, g in zip(np.floor(50+2*x), np.floor(30+2*y))]
    source = ColumnDataSource(data=dict(x=x,y=y,radius=radii))

    def makeplot():
        p = Figure(plot_width=800, plot_height=200,tools="hover,pan",title=None)
        p.scatter(x, y, radius=radii,
               fill_color=colors, fill_alpha=0.6,
               line_color='gray', source=source)
        hover = p.select(dict(type=HoverTool))
        hover.tooltips = OrderedDict([
            ("radius", "@radius")])
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p

    p1 = makeplot()
    p2 = makeplot()
    p3 = makeplot()
    p = GridPlot(children=[[p1],[p2],[p3]],title='test')
    js,html = embedPlot(p)
    return dict(figure=html)

@auth.requires_login()
def admin():
    """Settings"""
    parser,conffile = getConfig()
    options = dict(parser.items('base'))
    form = SQLFORM.dictform(options)
    if form.process().accepted:
        for i in dict(parser.items('base')):
            print i
            parser.set('base', i, form.vars[i])
        parser.write(open(conffile,'w'))
        response.flash='Saved'
        redirect(URL('default','admin'))
    return dict(form=form)

def about():
    msg = 'About this page'
    #fp = os.path.join(request.folder,'static/docs','about.txt')

    return dict(msg=msg)

def citation():
    return dict()

def help():
    msg = T('')
    return dict(msg=msg)

