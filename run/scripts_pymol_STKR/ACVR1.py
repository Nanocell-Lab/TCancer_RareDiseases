# Script de PyMOL para el gen ACVR1
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('R206H', 'resi 206')
cmd.show('sticks', 'R206H')
cmd.label('''byca(R206H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 206')

cmd.select('G356D', 'resi 356')
cmd.show('sticks', 'G356D')
cmd.label('''byca(G356D)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 356')

