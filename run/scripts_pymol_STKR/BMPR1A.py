# Script de PyMOL para el gen BMPR1A
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('R361*', 'resi 361')
cmd.show('sticks', 'R361*')
cmd.color('lightmagenta', 'resi 361')

cmd.select('R494*', 'resi 494')
cmd.show('sticks', 'R494*')
cmd.color('lightmagenta', 'resi 494')

cmd.select('R228*', 'resi 228')
cmd.show('sticks', 'R228*')
cmd.color('lightmagenta', 'resi 228')

cmd.select('W504*', 'resi 504')
cmd.show('sticks', 'W504*')
cmd.color('lightmagenta', 'resi 504')

cmd.select('R273*', 'resi 273')
cmd.show('sticks', 'R273*')
cmd.color('lightmagenta', 'resi 273')

cmd.select('Q22*', 'resi 22')
cmd.show('sticks', 'Q22*')
cmd.color('lightmagenta', 'resi 22')

