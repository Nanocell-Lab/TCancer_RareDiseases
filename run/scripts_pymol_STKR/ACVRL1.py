# Script de PyMOL para el gen ACVRL1
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('F425I', 'resi 425')
cmd.show('sticks', 'F425I')
cmd.color('lightmagenta', 'resi 425')

cmd.select('E224K', 'resi 224')
cmd.show('sticks', 'E224K')
cmd.color('lightmagenta', 'resi 224')

cmd.select('S462*', 'resi 462')
cmd.show('sticks', 'S462*')
cmd.color('lightmagenta', 'resi 462')

cmd.select('D330Y', 'resi 330')
cmd.show('sticks', 'D330Y')
cmd.color('lightmagenta', 'resi 330')

cmd.select('R411Q', 'resi 411')
cmd.show('sticks', 'R411Q')
cmd.color('lightmagenta', 'resi 411')

cmd.select('R200W', 'resi 200')
cmd.show('sticks', 'R200W')
cmd.color('lightmagenta', 'resi 200')

cmd.select('T277M', 'resi 277')
cmd.show('sticks', 'T277M')
cmd.color('lightmagenta', 'resi 277')

cmd.select('W472*', 'resi 472')
cmd.show('sticks', 'W472*')
cmd.color('lightmagenta', 'resi 472')

cmd.select('V441M', 'resi 441')
cmd.show('sticks', 'V441M')
cmd.color('lightmagenta', 'resi 441')

cmd.select('R374W', 'resi 374')
cmd.show('sticks', 'R374W')
cmd.color('lightmagenta', 'resi 374')

cmd.select('R484Q', 'resi 484')
cmd.show('sticks', 'R484Q')
cmd.color('lightmagenta', 'resi 484')

cmd.select('G350S', 'resi 350')
cmd.show('sticks', 'G350S')
cmd.color('lightmagenta', 'resi 350')

cmd.select('R386H', 'resi 386')
cmd.show('sticks', 'R386H')
cmd.color('lightmagenta', 'resi 386')

cmd.select('E236K', 'resi 236')
cmd.show('sticks', 'E236K')
cmd.color('lightmagenta', 'resi 236')

cmd.select('F425L', 'resi 425')
cmd.show('sticks', 'F425L')
cmd.color('lightmagenta', 'resi 425')

cmd.select('P378H', 'resi 378')
cmd.show('sticks', 'P378H')
cmd.color('lightmagenta', 'resi 378')

cmd.select('R144*', 'resi 144')
cmd.show('sticks', 'R144*')
cmd.color('lightmagenta', 'resi 144')

cmd.select('L381P', 'resi 381')
cmd.show('sticks', 'L381P')
cmd.color('lightmagenta', 'resi 381')

cmd.select('R67Q', 'resi 67')
cmd.show('sticks', 'R67Q')
cmd.color('lightmagenta', 'resi 67')

cmd.select('W420*', 'resi 420')
cmd.show('sticks', 'W420*')
cmd.color('lightmagenta', 'resi 420')

cmd.select('W406*', 'resi 406')
cmd.show('sticks', 'W406*')
cmd.color('lightmagenta', 'resi 406')

cmd.select('C443*', 'resi 443')
cmd.show('sticks', 'C443*')
cmd.color('lightmagenta', 'resi 443')

