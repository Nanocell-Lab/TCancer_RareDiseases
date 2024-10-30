# Script de PyMOL para el gen BMPR2
from pymol import cmd

cmd.bg_color(color='white')
cmd.orient()
cmd.set('cartoon_transparency', 0.4)
cmd.color('aquamarine', 'all')
cmd.select('R321*', 'resi 321')
cmd.show('sticks', 'R321*')
cmd.label('''byca(R321*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 321')

cmd.select('R213*', 'resi 213')
cmd.show('sticks', 'R213*')
cmd.label('''byca(R213*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 213')

cmd.select('S779*', 'resi 779')
cmd.show('sticks', 'S779*')
cmd.label('''byca(S779*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 779')

cmd.select('C420R', 'resi 420')
cmd.show('sticks', 'C420R')
cmd.label('''byca(C420R)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 420')

cmd.select('Y113C', 'resi 113')
cmd.show('sticks', 'Y113C')
cmd.label('''byca(Y113C)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 113')

cmd.select('R332*', 'resi 332')
cmd.show('sticks', 'R332*')
cmd.label('''byca(R332*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 332')

cmd.select('Y214*', 'resi 214')
cmd.show('sticks', 'Y214*')
cmd.label('''byca(Y214*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 214')

cmd.select('R147*', 'resi 147')
cmd.show('sticks', 'R147*')
cmd.label('''byca(R147*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 147')

cmd.select('C420Y', 'resi 420')
cmd.show('sticks', 'C420Y')
cmd.label('''byca(C420Y)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 420')

cmd.select('R491W', 'resi 491')
cmd.show('sticks', 'R491W')
cmd.label('''byca(R491W)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 491')

cmd.select('Q6*', 'resi 6')
cmd.show('sticks', 'Q6*')
cmd.label('''byca(Q6*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 6')

cmd.select('W16*', 'resi 16')
cmd.show('sticks', 'W16*')
cmd.label('''byca(W16*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 16')

cmd.select('W298*', 'resi 298')
cmd.show('sticks', 'W298*')
cmd.label('''byca(W298*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 298')

cmd.select('Y282*', 'resi 282')
cmd.show('sticks', 'Y282*')
cmd.label('''byca(Y282*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 282')

cmd.select('R211*', 'resi 211')
cmd.show('sticks', 'R211*')
cmd.label('''byca(R211*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 211')

cmd.select('N126S', 'resi 126')
cmd.show('sticks', 'N126S')
cmd.label('''byca(N126S)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 126')

cmd.select('Q450*', 'resi 450')
cmd.show('sticks', 'Q450*')
cmd.label('''byca(Q450*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 450')

cmd.select('R873*', 'resi 873')
cmd.show('sticks', 'R873*')
cmd.label('''byca(R873*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 873')

cmd.select('E98*', 'resi 98')
cmd.show('sticks', 'E98*')
cmd.label('''byca(E98*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 98')

cmd.select('W9*', 'resi 9')
cmd.show('sticks', 'W9*')
cmd.label('''byca(W9*)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 9')

cmd.select('E386G', 'resi 386')
cmd.show('sticks', 'E386G')
cmd.label('''byca(E386G)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 386')

cmd.select('Q82H', 'resi 82')
cmd.show('sticks', 'Q82H')
cmd.label('''byca(Q82H)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 82')

cmd.select('A24E', 'resi 24')
cmd.show('sticks', 'A24E')
cmd.label('''byca(A24E)''', 'oneletter+resi')
cmd.color('lightmagenta', 'resi 24')

