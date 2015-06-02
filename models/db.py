# -*- coding: utf-8 -*-

#########################################################################
## This scaffolding model makes your app work on Google App Engine too
## File is released under public domain and you can use without limitations
#########################################################################

## if SSL/HTTPS is properly configured and you want all HTTP requests to
## be redirected to HTTPS, uncomment the line below:
# request.requires_https()

if not request.env.web2py_runtime_gae:
    ## if NOT running on Google App Engine use SQLite or other DB
    db = DAL('sqlite://storage.sqlite')#,lazy_tables=True)
else:
    ## connect to Google BigTable (optional 'google:datastore://namespace')
    db = DAL('google:datastore')
    ## store sessions and tickets there
    session.connect(request, response, db = db)


## by default give a view/generic.extension to all actions from localhost
## none otherwise. a pattern can be 'controller/function.extension'
response.generic_patterns = ['*'] if request.is_local else []

#########################################################################
## Here is sample code if you need for
## - email capabilities
## - authentication (registration, login, logout, ... )
## - authorization (role based authorization)
## - services (xml, csv, json, xmlrpc, jsonrpc, amf, rss)
## - old style crud actions
## (more options discussed in gluon/tools.py)
#########################################################################

from gluon.tools import Auth, Crud, Service, PluginManager, prettydate
auth = Auth(db, hmac_key=Auth.get_or_create_key())
crud, service, plugins = Crud(db), Service(), PluginManager()

## create all tables needed by auth if not custom tables
auth.define_tables()

## configure email
mail=auth.settings.mailer
mail.settings.server = 'logging' or 'smtp.gmail.com:587'
mail.settings.sender = 'you@gmail.com'
mail.settings.login = 'username:password'

## configure auth policy
auth.settings.registration_requires_verification = False
auth.settings.registration_requires_approval = False
auth.settings.reset_password_requires_verification = True

#auth.settings.register_onaccept.append(lambda form:
#    mail.send(to='ad...@emai.com',subject='New user registered for %s application' % (request.application),
#                message="new user email is %s" % (form.vars.email)))

## Define your tables below (or better in another model file) for example

db.define_table('genomes',
        Field('name','string',unique=True,requires=IS_NOT_EMPTY()),
        Field('filename'),
        Field('format','string'),
		Field('file','upload'))#,migrate=False)

db.define_table('predictions',
        Field('identifier','string',required=True),
        Field('user','string',requires=IS_NOT_EMPTY()),
        Field('description','string',requires=IS_NOT_EMPTY()))
		#Field('path','string',unique=True,requires=IS_NOT_EMPTY()))

db.define_table('proteins',
        Field('genome','string'),
		Field('locus_tag','string'),
        Field('type','string'),
		Field('protein_id','string'),
        Field('product','string'),
		Field('gene','string'),
		Field('start','integer'),
        Field('end','integer'),
        Field('length','integer'),
		Field('db_xref','string'),
		Field('note','text'),
		Field('pseudo','string'),
		Field('translation','text')	)

from gluon.custom_import import track_changes; track_changes(True)
import sys, os
import numpy as np
import pylab as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import cStringIO


