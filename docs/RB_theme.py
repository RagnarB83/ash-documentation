# -*- coding: utf-8 -*-
"""
    pygments.styles.RBthemeStyle
    ~~~~~~~~~~~~~~~~~~~~~~~
    RB. MOdification to Monokai to mimic okaidia
    Mimic the Monokai color scheme. Based on tango.py.

    http://www.monokai.nl/blog/2006/07/15/textmate-color-theme/

    :copyright: Copyright 2006-2019 by the Pygments team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

#RB. Not

from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, Text, \
     Number, Operator, Generic, Whitespace, Punctuation, Other, Literal

class Rb_ThemeStyle(Style):
    """
    This style mimics the Monokai color scheme.
    """

    #background_color = "#272822"
    background_color = "#272822"
    highlight_color = "#49483e"

    styles = {
        # No corresponding class for the following:
        Text:                      "#f8f8f2", # class:  ''
        Whitespace:                "",        # class: 'w'
        Error:                     "#960050 bg:#1e0010", # class: 'err'
        Other:                     "",        # class 'x'

        Comment:                   "#708090", # class: 'c'
        Comment.Multiline:         "",        # class: 'cm'
        Comment.Preproc:           "",        # class: 'cp'
        Comment.Single:            "",        # class: 'c1'
        Comment.Special:           "",        # class: 'cs'

        Keyword:                   "#66d9ef", # class: 'k'
        Keyword.Constant:          "",        # class: 'kc'
        Keyword.Declaration:       "",        # class: 'kd'
        #Keyword.Namespace:         "#f92672", # class: 'kn'
        Keyword.Namespace:         "#66d9ef", # class: 'kn'
        Keyword.Pseudo:            "",        # class: 'kp'
        Keyword.Reserved:          "",        # class: 'kr'
        Keyword.Type:              "",        # class: 'kt'

        #Operator:                  "#f92672", # class: 'o'
        Operator.Word:             "",        # class: 'ow' - like keywords

        #Punctuation:               "#f8f8f2", # class: 'p'

        Name:                      "#f8f8f2", # class: 'n'
        Name.Attribute:            "#a6e22e", # class: 'na' - to be revised
        Name.Builtin:              "",        # class: 'nb'
        Name.Builtin.Pseudo:       "",        # class: 'bp'
        Name.Class:                "#a6e22e", # class: 'nc' - to be revised
        Name.Constant:             "#66d9ef", # class: 'no' - to be revised
        Name.Decorator:            "#a6e22e", # class: 'nd' - to be revised
        Name.Entity:               "",        # class: 'ni'
        Name.Exception:            "#a6e22e", # class: 'ne'
        Name.Function:             "#a6e22e", # class: 'nf'
        Name.Property:             "",        # class: 'py'
        Name.Label:                "",        # class: 'nl'
        Name.Namespace:            "",        # class: 'nn' - to be revised
        Name.Other:                "#a6e22e", # class: 'nx'
        #Name.Tag:                  "#f92672", # class: 'nt' - like a keyword
        Name.Variable:             "",        # class: 'nv' - to be revised
        Name.Variable.Class:       "",        # class: 'vc' - to be revised
        Name.Variable.Global:      "",        # class: 'vg' - to be revised
        Name.Variable.Instance:    "",        # class: 'vi' - to be revised

        Number:                    "#ae81ff", # class: 'm'
        Number.Float:              "",        # class: 'mf'
        Number.Hex:                "",        # class: 'mh'
        Number.Integer:            "",        # class: 'mi'
        Number.Integer.Long:       "",        # class: 'il'
        Number.Oct:                "",        # class: 'mo'

        Literal:                   "#ae81ff", # class: 'l'
        Literal.Date:              "#e6db74", # class: 'ld'

        #String:                    "#e6db74", # class: 's'
        #String:                    "#e6db74", # class: 's'
        String:                    "#a6e22e", # class: 's'
        String.Backtick:           "",        # class: 'sb'
        String.Char:               "",        # class: 'sc'
        String.Doc:                "",        # class: 'sd' - like a comment
        String.Double:             "",        # class: 's2'
        String.Escape:             "#ae81ff", # class: 'se'
        String.Heredoc:            "",        # class: 'sh'
        String.Interpol:           "",        # class: 'si'
        String.Other:              "",        # class: 'sx'
        String.Regex:              "",        # class: 'sr'
        String.Single:             "",        # class: 's1'
        String.Symbol:             "",        # class: 'ss'


        Generic:                   "",        # class: 'g'
        #Generic.Deleted:           "#f92672", # class: 'gd',
        Generic.Emph:              "italic",  # class: 'ge'
        Generic.Error:             "",        # class: 'gr'
        Generic.Heading:           "",        # class: 'gh'
        Generic.Inserted:          "#a6e22e", # class: 'gi'
        Generic.Output:            "#66d9ef", # class: 'go'
        #Generic.Prompt:            "bold #f92672", # class: 'gp'
        Generic.Strong:            "bold",    # class: 'gs'
        Generic.Subheading:        "#75715e", # class: 'gu'
        Generic.Traceback:         "",        # class: 'gt'
    }
