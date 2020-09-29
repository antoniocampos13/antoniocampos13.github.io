#!/usr/bin/env python
# -*- coding: utf-8 -*- #

AUTHOR = 'Antonio Victor Campos Coelho'
SITENAME = "Antonio's Portfolio"
RELATIVE_URLS = False
SITEURL = 'https://antoniocampos13.github.io/'
SITETITLE = SITENAME
SITEDESCRIPTION = "Data Science Portfolio by Antonio Victor Campos Coelho"
SITESUBTITLE = 'PhD in Genetics'
SITELOGO = 'https://avatars.githubusercontent.com/antoniocampos13'
BROWSER_COLOR = '#333333'
PYGMENTS_STYLE = 'monokai'

USE_FOLDER_AS_CATEGORY = True
MAIN_MENU = True
HOME_HIDE_TAGS = False

DELETE_OUTPUT_DIRECTORY = False
PATH = 'content'
STATIC_PATHS = ['images','static','pages']

THEME = 'Flex'

TIMEZONE = 'America/Recife'

DEFAULT_LANG = 'en'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Blogroll
LINKS = (
('Brazilian Lattes CV', 'http://lattes.cnpq.br/2986394950644755'),
('ORCID', 'https://orcid.org/0000-0003-2143-9701'),
('Researcher ID/Publons', 'https://publons.com/researcher/2414076/antonio-victor-c-coelho/'),
)

# Social widget
SOCIAL = (
('github', 'https://github.com/antoniocampos13'),
('linkedin', 'https://www.linkedin.com/in/antonio-coelho-9aa338164/'),
)

DEFAULT_PAGINATION = 30

CC_LICENSE = {
    'name': 'Creative Commons Attribution-ShareAlike',
    'version': '4.0',
    'slug': 'by-sa'
}

COPYRIGHT_NAME = AUTHOR
COPYRIGHT_YEAR = 2020