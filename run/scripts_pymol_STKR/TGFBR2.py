# Script de PyMOL para el gen TGFBR2
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('D471N', 'resi 471')
cmd.show('sticks', 'D471N')
cmd.label('''byca(D471N)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 471')

cmd.select('D446N', 'resi 446')
cmd.show('sticks', 'D446N')
cmd.label('''byca(D446N)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 446')

cmd.select('R485C', 'resi 485')
cmd.show('sticks', 'R485C')
cmd.label('''byca(R485C)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 485')

cmd.select('R460C', 'resi 460')
cmd.show('sticks', 'R460C')
cmd.label('''byca(R460C)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 460')

cmd.select('R553C', 'resi 553')
cmd.show('sticks', 'R553C')
cmd.label('''byca(R553C)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 553')

cmd.select('R528C', 'resi 528')
cmd.show('sticks', 'R528C')
cmd.label('''byca(R528C)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 528')

cmd.select('R497*', 'resi 497')
cmd.show('sticks', 'R497*')
cmd.label('''byca(R497*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 497')

cmd.select('R522*', 'resi 522')
cmd.show('sticks', 'R522*')
cmd.label('''byca(R522*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 522')

cmd.select('R562C', 'resi 562')
cmd.show('sticks', 'R562C')
cmd.label('''byca(R562C)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 562')

cmd.select('R537C', 'resi 537')
cmd.show('sticks', 'R537C')
cmd.label('''byca(R537C)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 537')

cmd.select('R553H', 'resi 553')
cmd.show('sticks', 'R553H')
cmd.label('''byca(R553H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 553')

cmd.select('R528H', 'resi 528')
cmd.show('sticks', 'R528H')
cmd.label('''byca(R528H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 528')

cmd.select('E524*', 'resi 524')
cmd.show('sticks', 'E524*')
cmd.label('''byca(E524*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 524')

cmd.select('E499*', 'resi 499')
cmd.show('sticks', 'E499*')
cmd.label('''byca(E499*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 499')

cmd.select('W287R', 'resi 287')
cmd.show('sticks', 'W287R')
cmd.label('''byca(W287R)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 287')

cmd.select('W312R', 'resi 312')
cmd.show('sticks', 'W312R')
cmd.label('''byca(W312R)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 312')

cmd.select('K277N', 'resi 277')
cmd.show('sticks', 'K277N')
cmd.label('''byca(K277N)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 277')

cmd.select('K302N', 'resi 302')
cmd.show('sticks', 'K302N')
cmd.label('''byca(K302N)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 302')

cmd.select('R495*', 'resi 495')
cmd.show('sticks', 'R495*')
cmd.label('''byca(R495*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 495')

cmd.select('R520*', 'resi 520')
cmd.show('sticks', 'R520*')
cmd.label('''byca(R520*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 520')

cmd.select('R485H', 'resi 485')
cmd.show('sticks', 'R485H')
cmd.label('''byca(R485H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 485')

cmd.select('R460H', 'resi 460')
cmd.show('sticks', 'R460H')
cmd.label('''byca(R460H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 460')

cmd.select('R537H', 'resi 537')
cmd.show('sticks', 'R537H')
cmd.label('''byca(R537H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 537')

cmd.select('D549N', 'resi 549')
cmd.show('sticks', 'D549N')
cmd.label('''byca(D549N)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 549')

cmd.select('D524N', 'resi 524')
cmd.show('sticks', 'D524N')
cmd.label('''byca(D524N)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 524')

cmd.select('W546*', 'resi 546')
cmd.show('sticks', 'W546*')
cmd.label('''byca(W546*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 546')

cmd.select('W521*', 'resi 521')
cmd.show('sticks', 'W521*')
cmd.label('''byca(W521*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 521')

