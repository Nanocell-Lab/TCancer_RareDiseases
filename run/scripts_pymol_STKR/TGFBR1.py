# Script de PyMOL para el gen TGFBR1
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('R410Q', 'resi 410')
cmd.show('sticks', 'R410Q')
cmd.label('''byca(R410Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 410')

cmd.select('R487Q', 'resi 487')
cmd.show('sticks', 'R487Q')
cmd.label('''byca(R487Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 487')

cmd.select('R491Q', 'resi 491')
cmd.show('sticks', 'R491Q')
cmd.label('''byca(R491Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 491')

cmd.select('S164L', 'resi 164')
cmd.show('sticks', 'S164L')
cmd.label('''byca(S164L)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 164')

cmd.select('S241L', 'resi 241')
cmd.show('sticks', 'S241L')
cmd.label('''byca(S241L)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 241')

cmd.select('S245L', 'resi 245')
cmd.show('sticks', 'S245L')
cmd.label('''byca(S245L)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 245')

cmd.select('R410W', 'resi 410')
cmd.show('sticks', 'R410W')
cmd.label('''byca(R410W)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 410')

cmd.select('R487W', 'resi 487')
cmd.show('sticks', 'R487W')
cmd.label('''byca(R487W)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 487')

cmd.select('R491W', 'resi 491')
cmd.show('sticks', 'R491W')
cmd.label('''byca(R491W)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 491')

cmd.select('W200R', 'resi 200')
cmd.show('sticks', 'W200R')
cmd.label('''byca(W200R)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 200')

cmd.select('W277R', 'resi 277')
cmd.show('sticks', 'W277R')
cmd.label('''byca(W277R)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 277')

cmd.select('W281R', 'resi 281')
cmd.show('sticks', 'W281R')
cmd.label('''byca(W281R)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 281')

cmd.select('R157*', 'resi 157')
cmd.show('sticks', 'R157*')
cmd.label('''byca(R157*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 157')

cmd.select('R161*', 'resi 161')
cmd.show('sticks', 'R161*')
cmd.label('''byca(R161*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 161')

cmd.select('L277P', 'resi 277')
cmd.show('sticks', 'L277P')
cmd.label('''byca(L277P)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 277')

cmd.select('L354P', 'resi 354')
cmd.show('sticks', 'L354P')
cmd.label('''byca(L354P)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 354')

cmd.select('L358P', 'resi 358')
cmd.show('sticks', 'L358P')
cmd.label('''byca(L358P)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 358')

