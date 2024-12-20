# Script de PyMOL para el gen BMPR1A
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('Q22*', 'resi 22')
cmd.show('sticks', 'Q22*')
cmd.label('''byca(Q22*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 22')

cmd.select('R228*', 'resi 228')
cmd.show('sticks', 'R228*')
cmd.label('''byca(R228*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 228')

cmd.select('R273*', 'resi 273')
cmd.show('sticks', 'R273*')
cmd.label('''byca(R273*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 273')

cmd.select('W504*', 'resi 504')
cmd.show('sticks', 'W504*')
cmd.label('''byca(W504*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 504')

cmd.select('R361*', 'resi 361')
cmd.show('sticks', 'R361*')
cmd.label('''byca(R361*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 361')

cmd.select('R494*', 'resi 494')
cmd.show('sticks', 'R494*')
cmd.label('''byca(R494*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 494')

