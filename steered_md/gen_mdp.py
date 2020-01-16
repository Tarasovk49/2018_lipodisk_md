
from sys import argv
from getopt import getopt

opts, args = getopt(argv[1:], 'm:i:o:r:e:')

input_mdp = 'config/lipodisk.mdp'
index = 'index.ndx'
output_mdp = 'config/lipodisk_flatbot.mdp'
R = 4
E = 50

for o, a in opts:
    if o == '-m':
        # input mdp file (everything before pull options)
        input_mdp = a
    if o == '-i':
        # index file with groups POL1, POL2, .., MEMBRANE
        # POLs will be pulled towards MEMBRANE
        index = a
    if o == '-o':
        # output mdp file with pull options
        output_mdp = a
    if o == '-r':
        # Forces are applied while distances are higher than R nm
        R = int(a)
    if o == '-e':
        # Force constant, kJ mol^-1 nm^-2
        E = int(a)

def gen_mdp(input_mdp, index, output_mdp):
    """
    Adds pull options to mdp file. Center of mass of each polymer molecule is pulled
    towards center of mass of membrane by applying harmonic potential until distance become @R nm
    """
    polymer_mol_count=0
    with open(index) as ndx:
        for line in ndx.readlines():
            if line[0:5] == '[ POL':
                polymer_mol_count+=1
    with open(output_mdp,'w') as mdp:
        with open(input_mdp) as main_mdp:
            for string in main_mdp.readlines():
                mdp.write(string)
                
            mdp.write('; Pull code'+'\n')
            mdp.write('pull                    = yes'+'\n')
            mdp.write('pull_ngroups            = {0}'.format(polymer_mol_count+1)+'\n')
            mdp.write('pull_ncoords            = {0}'.format(polymer_mol_count)+'\n')
            mdp.write('pull_group1_name        = MEMBRANE'+'\n')
            
            for i in range(1,1+polymer_mol_count):
                mdp.write('pull_group{0}_name        = POL{1}'.format((i+1),i)+'\n')
                mdp.write('pull_coord{0}_type        = flat-bottom   ; harmonic biasing force pulling only when distance > pull_coord_init'.format(i)+'\n')
                mdp.write('pull_coord{0}_geometry    = distance      ; '.format(i)+'\n')
                mdp.write('pull_coord{0}_groups      = 1 {1}'.format(i,(i+1))+'\n')
                mdp.write('pull_coord{0}_dim         = Y Y N'.format(i)+'\n')
                mdp.write('pull_coord{0}_init        = {1}           ; distance decreases while it is higher than {1} nm'.format(i,R)+'\n')
                mdp.write('pull_coord{0}_k           = {1}            ; kJ mol^-1 nm^-2'.format(i,E)+'\n')
                mdp.write('\n')
                                
gen_mdp(input_mdp, index, output_mdp)
