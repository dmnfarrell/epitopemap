# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## Customize your APP title, subtitle and menus here
#########################################################################

response.title = ' '.join(word.capitalize() for word in request.application.split('_'))
response.subtitle = T('Visualise MHC binding predictions')

## read more at http://dev.w3.org/html5/markup/meta.name.html
response.meta.author = 'Damien Farrell <damien.farrell@ucd.ie>'
response.meta.description = 'Epitopemap application'
response.meta.generator = 'Web2py Web Framework'
response.meta.copyright = 'Copyright 2014-'

## your http://google.com/analytics id
response.google_analytics_id = None

#########################################################################
## this is the main application menu add/remove items as required
#########################################################################

response.menu = [
    (T('Home'), False, URL('default','index'), []),
    (T('Help'), False, URL('default','help'), []),
    (T('About'), False, URL('default','about'), [])
    ]

DEVELOPMENT_MENU = True

#########################################################################
## provide shortcuts for development. remove in production
#########################################################################

def _():
    # shortcuts
    app = request.application
    ctr = request.controller
    # useful links to internal and external resources
    response.menu+=[
        (SPAN('web2py',_style='color:yellow'),False, None, [
                (T('This App'),False,URL('admin','default','design/%s' % app), [
                        (T('Controller'),False,
                         URL('admin','default','edit/%s/controllers/%s.py' % (app,ctr))),
                        (T('View'),False,
                         URL('admin','default','edit/%s/views/%s' % (app,response.view))),
                        (T('Layout'),False,
                         URL('admin','default','edit/%s/views/layout.html' % app)),
                        (T('Stylesheet'),False,
                         URL('admin','default','edit/%s/static/css/web2py.css' % app)),
                        (T('DB Model'),False,
                         URL('admin','default','edit/%s/models/db.py' % app)),
                        (T('Menu Model'),False,
                         URL('admin','default','edit/%s/models/menu.py' % app)),
                        (T('Database'),False, URL(app,'appadmin','index')),
                        (T('Errors'),False, URL('admin','default','errors/' + app)),
                        (T('About'),False, URL('admin','default','about/' + app)),
                        ])
                ]
         )]
_()

if "auth" in locals(): auth.wikimenu()
