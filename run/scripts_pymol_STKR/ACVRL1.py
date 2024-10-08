# Script de PyMOL para el gen ACVRL1
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('F425I', 'resi 425')
cmd.show('sticks', 'F425I')
cmd.label('''byca(F425I)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 425')

cmd.select('E224K', 'resi 224')
cmd.show('sticks', 'E224K')
cmd.label('''byca(E224K)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 224')

cmd.select('S462*', 'resi 462')
cmd.show('sticks', 'S462*')
cmd.label('''byca(S462*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 462')

cmd.select('D330Y', 'resi 330')
cmd.show('sticks', 'D330Y')
cmd.label('''byca(D330Y)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 330')

cmd.select('R411Q', 'resi 411')
cmd.show('sticks', 'R411Q')
cmd.label('''byca(R411Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 411')

cmd.select('R200W', 'resi 200')
cmd.show('sticks', 'R200W')
cmd.label('''byca(R200W)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 200')

cmd.select('T277M', 'resi 277')
cmd.show('sticks', 'T277M')
cmd.label('''byca(T277M)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 277')

cmd.select('W472*', 'resi 472')
cmd.show('sticks', 'W472*')
cmd.label('''byca(W472*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 472')

cmd.select('V441M', 'resi 441')
cmd.show('sticks', 'V441M')
cmd.label('''byca(V441M)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 441')

cmd.select('R374W', 'resi 374')
cmd.show('sticks', 'R374W')
cmd.label('''byca(R374W)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 374')

cmd.select('R484Q', 'resi 484')
cmd.show('sticks', 'R484Q')
cmd.label('''byca(R484Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 484')

cmd.select('G350S', 'resi 350')
cmd.show('sticks', 'G350S')
cmd.label('''byca(G350S)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 350')

cmd.select('R386H', 'resi 386')
cmd.show('sticks', 'R386H')
cmd.label('''byca(R386H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 386')

cmd.select('E236K', 'resi 236')
cmd.show('sticks', 'E236K')
cmd.label('''byca(E236K)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 236')

cmd.select('F425L', 'resi 425')
cmd.show('sticks', 'F425L')
cmd.label('''byca(F425L)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 425')

cmd.select('P378H', 'resi 378')
cmd.show('sticks', 'P378H')
cmd.label('''byca(P378H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 378')

cmd.select('R144*', 'resi 144')
cmd.show('sticks', 'R144*')
cmd.label('''byca(R144*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 144')

cmd.select('L381P', 'resi 381')
cmd.show('sticks', 'L381P')
cmd.label('''byca(L381P)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 381')

cmd.select('R67Q', 'resi 67')
cmd.show('sticks', 'R67Q')
cmd.label('''byca(R67Q)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 67')

cmd.select('W406*', 'resi 406')
cmd.show('sticks', 'W406*')
cmd.label('''byca(W406*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 406')

cmd.select('W420*', 'resi 420')
cmd.show('sticks', 'W420*')
cmd.label('''byca(W420*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 420')

cmd.select('C443*', 'resi 443')
cmd.show('sticks', 'C443*')
cmd.label('''byca(C443*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 443')

