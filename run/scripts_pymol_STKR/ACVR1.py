# Script de PyMOL para el gen ACVR1
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('R206H', 'resi 206')
cmd.wizard('mutagenesis')
cmd.get_wizard().set_mode('HIS')
cmd.get_wizard().do_select('resi 206')
cmd.get_wizard().apply()
cmd.show('sticks', 'resi 206')
cmd.select('HIS206', 'resi 206')
cmd.label('''byca(HIS)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 206')

cmd.select('G356D', 'resi 356')
cmd.wizard('mutagenesis')
cmd.get_wizard().set_mode('ASP')
cmd.get_wizard().do_select('resi 356')
cmd.get_wizard().apply()
cmd.show('sticks', 'resi 356')
cmd.select('ASP356', 'resi 356')
cmd.label('''byca(ASP)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 356')

