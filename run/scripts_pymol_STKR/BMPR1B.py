# Script de PyMOL para el gen BMPR1B
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('R486Q', 'resi 486')
cmd.show('sticks', 'R486Q')
cmd.label('''byca(R486Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 486')

cmd.select('R516Q', 'resi 516')
cmd.show('sticks', 'R516Q')
cmd.label('''byca(R516Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 516')

