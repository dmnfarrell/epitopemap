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
response.meta.version = '1.0'

## your http://google.com/analytics id
response.google_analytics_id = None

#########################################################################
## this is the main application menu add/remove items as required
#########################################################################

response.menu = [
    (T('Home'), False, URL('default','index'), []),
    (T('Help'), False, URL('default','help'), []),
    ]

DEVELOPMENT_MENU = True

#########################################################################
## provide shortcuts for development. remove in production
#########################################################################

def _():
    # shortcuts
    app = request.application
    ctr = request.controller
    response.menu+=[
        (SPAN('Menu',_style='color:yellow'),False, None, [
            (T('View predictions'),False,URL('default','quickview')),
            (T('Search'),False,URL('default','search')),
            (T('Submit job'),False,URL('default','submit')),
            (T('Configuration'),False,URL('default','admin')),
            (T('About'),False,URL('default','about')),
            (T('Site admin'),False,URL('admin','default','design/%s' % app)),]
         )]
_()

if "auth" in locals(): auth.wikimenu()
