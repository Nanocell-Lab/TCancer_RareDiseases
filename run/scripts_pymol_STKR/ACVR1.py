# Script de PyMOL para el gen ACVR1
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('G328V', 'resi 328')
cmd.show('sticks', 'G328V')
cmd.color('lightmagenta', 'resi 328')

cmd.select('R206H', 'resi 206')
cmd.show('sticks', 'R206H')
cmd.color('lightmagenta', 'resi 206')

cmd.select('G356D', 'resi 356')
cmd.show('sticks', 'G356D')
cmd.color('lightmagenta', 'resi 356')

