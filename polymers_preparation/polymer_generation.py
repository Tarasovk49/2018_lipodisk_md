from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import copolymer
import numpy as np
import shutil, os, sys
from getopt import getopt

##########################################################################################################################
#####################          Script builds styrol-diisobutylene-maleic acid copolymers.            #####################
##########################################################################################################################
#                                                                                                                        #
# The options provided are (with defaults):                                                                              #
#                                                                                                                        #
# '-n': (100) number of polymers to be generated. Is taken into account only with '--random'='yes'                       #
#                                                                                                                        #
# '--length': (36.0) mean length of polymer molecule in monomers. Is taken into account only with '--random'='yes'       #
#                                                                                                                        #
# '--rmsd': (3.0) RMSD of lengths of polymer molecules in monomers. Is taken into account only with '--random'='yes'     #
#                                                                                                                        #
# '--protstate': (3) Defines the charge on MA monomers (acidity of media).                                               #
#                                                                                                                        #
#               Protonation states of maleic acid residues:                                                              #
#               protstate    1     2     3     4     5     6     7                                                       #
#               pH           5     6     7     8     9    10   >10                                                       #
#               charge    -0.3  -0.5  -1.0  -1.2  -1.7  -1.9  -2.0                                                       #
#                                                                                                                        #
# '--pm': (0.25)              occurence of maleic acid monomer. Is taken into account only with '--random'='yes'         #
#                                                                                                                        #
# '--ps': (0.75)              occurence of styrol monomer. Is taken into account only with '--random'='yes'              #
#                                                                                                                        #
# '--pd': (0.0)              occurence of diisobutylene monomer. Is taken into account only with '--random'='yes'        #
#                                                                                                                        #
# '--random': (yes) You may support a text file with sequences ('--random'='no'), write them manually ('--random'='no'), #
#                   or specify all the flags listed above to generate random polymers ('--random'='yes')                 #
#                                                                                                                        #
# '--seq_file': (not_specified) You may support a txt file with sequences of polymer molecules.                          #
#                               Every string is the sequence of one polymer moleculeIt must look like:                   #
#                                                                                                                        #
#                         MSMSMMMMMMSSSSSSMSMS                                                                           #
#                         SSMSSMSSSMSMSMSSSSMSSM                                                                         #
#                         SSMMSMSMSMSMSSMMS                                                                              #
#                         ...                                                                                            #
#                                                                                                                        #
# '--dir_name': (polymer_molecules) name of folder to store pdb files for polymer molecules                              #
#                                                                                                                        #
# '--out_filename': (polymer) prefix of name of pdb files for polymer molecules                                          #
#                                                                                                                        #
#                                                                                                                        #
# @position variable which is passed to monomer build functions declares                                                 #
# one of three positions monomer can reside: first(@position=1), last(@position=2) and middle(@position=0).              #
# The difference between positions is that in case of first and last monomers terminal                                   #
# hydrogens must not be removed during "polimerization".                                                                 #
# There are 8 types of monomers:                                                                                         #
# 4 = (2 different states of protonation) x (2 different carboxyls of maleic acid),                                      #
# 2 for diisobutylene monomers with radical group oriented closer to head or tail of polymer molecule,                   #
# 2 for styrene monomers with benzene group oriented closer to head or tail of polymer molecule.                         #
##########################################################################################################################


def dib_monomer(position=0):
    try:
        s = system.read_mol('topology_generalized/DIB1.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    
    c1 = s.particles[3]
    c2 = s.particles[7]
    c1.linker = 'head'
    c2.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'DB1'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c2.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c2 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()

    return s

def dib_monomer_reverted(position=0):
    try:
        s = system.read_mol('topology_generalized/DIB2.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    
    c1 = s.particles[7]
    c2 = s.particles[3]
    c1.linker = 'head'
    c2.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'DB2'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c2.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c2 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()

    return s
        
def st_monomer(position=0):
    try:
        s = system.read_mol('topology_generalized/STY1.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    

    c1 = s.particles[1]
    c5 = s.particles[5]
    c1.linker = 'head'
    c5.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'ST1'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c5.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c5 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()
 
    return s

def st_monomer_reverted(position=0):
    try:
        s = system.read_mol('topology_generalized/STY2.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    

    c5 = s.particles[1]
    c1 = s.particles[5]
    c1.linker = 'head'
    c5.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'ST2'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c5.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c5 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()
 
    return s

def ma2_monomer(position=0):
    try:
        s = system.read_mol('topology_generalized/MA2.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    
    c1 = s.particles[5]
    c2 = s.particles[6]
    c1.linker = 'head'
    c2.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'MA2'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c2.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c2 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()

    return s

def mal_monomer(position=0):
    try:
        s = system.read_mol('topology_generalized/MAL.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    
    c1 = s.particles[5]
    c2 = s.particles[6]
    c1.linker = 'head'
    c2.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'MAL'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c2.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c2 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()

    return s
    
def mar_monomer(position=0):
    try:
        s = system.read_mol('topology_generalized/MAR.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    
    c1 = s.particles[5]
    c2 = s.particles[6]
    c1.linker = 'head'
    c2.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'MAR'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c2.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c2 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()

    return s

def mad_monomer(position=0):
    try:
        s = system.read_mol('topology_generalized/MAD.mol')
    except:
        pass
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    
    c1 = s.particles[5]
    c2 = s.particles[6]
    c1.linker = 'head'
    c2.linker = 'tail'
    
    for p in s.particles:
        p.resid = 0
        p.resname = 'MAD'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            if position != 2:
                s.particles.remove(pb.tag, update=False)
            break
        
    for b in c2.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c2 else b.b
            if position != 1:
                s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')   
    s.add_particle_bonding()

    return s



if __name__ == '__main__':  
    
    monomers={'db' : dib_monomer(),
             'dbf' : dib_monomer(1),
             'dbl' : dib_monomer(2),
             'dbre' : dib_monomer_reverted(),
             'dbrf' : dib_monomer_reverted(1),
             'dbrl' : dib_monomer_reverted(2),
             'ma2' : ma2_monomer(),
             'ma2f' : ma2_monomer(1),
             'ma2l' : ma2_monomer(2),
             'mal' : mal_monomer(),
             'malf' : mal_monomer(1),
             'mall' : mal_monomer(2),
             'mar' : mar_monomer(),
             'marf' : mar_monomer(1),
             'marl' : mar_monomer(2),
             'mad' : mad_monomer(),
             'madf' : mad_monomer(1),
             'madl' : mad_monomer(2),
             'st' : st_monomer(),
             'stf' : st_monomer(1),
             'stl' : st_monomer(2),
             'stre' : st_monomer_reverted(),
             'strf' : st_monomer_reverted(1),
             'strl' : st_monomer_reverted(2)}
            

        
    opts, args = getopt(sys.argv[1:], 'n:',longopts=['length=','protstate=','rmsd=','out_filename=','pm=','ps=','pd=','random=','seq_file=','dir_name='])
    
    ####################################################
    ##################### DEFAULTS #####################
    ####################################################
    
    # Option to specify if you want to build random polymers or polymers of particular sequences
    # supplied from txt file or manually
    random = 'yes'
    
    
    # You may want to build polymers with exact set of sequences. Then provide txt file
    # and support it to '--seq_file' option or here. It must look like:
    # '''
    # MSMSMMMMMMSSSSSSMSMS
    # SSMSSMSSSMSMSMSSSSMSSM
    # SSMMSMSMSMSMSSMMS
    # ...
    # '''
    # Every string is the sequence of one polymer molecule
    seq_file='not_specified'
    
    
    # Number of polymer molecules to be generated, mean of their length, rmsd of length,
    # protonation state of MA monomers
    n = 100
    meanlen = 36.0
    rmsd = 3.0
    protstate = 3

    # These are the occurrencies of monomers (pm - maleic acid, pd - diisobutylene, ps - styrol)
    pm = 0.25
    pd = 0.0
    ps = 0.75
    
    # prefix of name of pdb files for polymer molecules and name of folder to store them
    out_filename = 'polymer'
    dir_name = 'polymer_molecules'
    
    ####################################################
    ################# END OF DEFAULTS ##################
    ####################################################
    
    # Default values are substituted here if options were specified    
    for o, a in opts:
        if o == '-n':
            n = int(a)
        if o == '--length':
            meanlen = float(a)
        if o == '--rmsd':
            rmsd = float(a)
        if o == '--protstate':
            protstate = int(a)
        if o == '--pm':
            pm = float(a)
        if o == '--ps':
            ps = float(a)
        if o == '--pd':
            pd = float(a)
        if o == '--random':
            random = str(a)
        if o == '--seq_file':
            seq_file = str(a)
        if o == '--dir_name':
            dir_name = str(a)
        if o == '--out_filename':
            out_filename = str(a)
    
    
    
    # Create folder for polymer molecules to store
    if os.path.isdir(dir_name):
        shutil.rmtree(dir_name)
        os.mkdir(dir_name)
    else:
        os.mkdir(dir_name)
    
    # During polymerization great amount of data is generated. It makes the main folder
    # look like a Trash. So we build polymer inside tmp folder and then remove it.
    if os.path.isdir('tmp'):
        shutil.rmtree('tmp')
        os.mkdir('tmp')
        os.chdir('tmp')
    else:
        os.mkdir('tmp')
        os.chdir('tmp')
    

    
    # Occurencies of monomers
    print "Normalizing occurrencies..."
    s = pm + ps + pd
    pm /= s
    ps /= s
    pd /= s
    print "Normalized occurrencies are:"
    print('ps = '+str(ps)+'; pm = '+str(pm)+'; pd = '+str(pd))


    db = monomers['db']
    dbf = monomers['dbf']
    dbl = monomers['dbl']
    
    dbre = monomers['dbre']
    dbrf = monomers['dbrf']
    dbrl = monomers['dbrl']

    ma2 = monomers['ma2']
    ma2f = monomers['ma2f']
    ma2l = monomers['ma2l']
    
    mal = monomers['mal']
    malf = monomers['malf']
    mall = monomers['mall']
    
    mar = monomers['mar']
    marf = monomers['marf']
    marl = monomers['marl']
    
    mad = monomers['mad']
    madf = monomers['madf']
    madl = monomers['madl']    
    
    st = monomers['st']
    stf = monomers['stf']   
    stl = monomers['stl']
    
    stre = monomers['stre']
    strf = monomers['strf']
    strl = monomers['strl']   
    f = forcefield.Dreiding()
    
    
    # This part helps you not to make mistake
    
    if (random=='yes') and (seq_file != 'not_specified'):
        print('\nSeems like you have just put it up to eleven!!')
        print('\nYou specified file with sequences '+seq_file+' and nonetheless want to build random polymers!\n')
        from_file='thisstringdoesntmeananything'
        while True:
            if not ((from_file == 'yes') or (from_file == 'no')):
                from_file = raw_input('Do you want to build polymer molecules with sequences from '+seq_file+'? Type "yes" or "no". \n')
                continue
            else:
                break
        if from_file=='yes':
            random='no'
        elif from_file=='no':
            seq_file ='not_specified'
            
  
    # NOTE that from this point 'ps' and 'pd' are not the probabilities to have styrol or diisobutylene
    # in any position, but probabilities to have left-oriented or right-oriented styrol\diisobutylene in that position!!
    ps = ps/2
    pd = pd/2  
  
            
            
    if seq_file != 'not_specified':
        ###############
        ## FROM FILE ##
        ###############
        with open('../'+seq_file) as pattern:
            for molecule_number, string in enumerate(pattern):
                string = string.strip('\n')
                length = len(string)
                seq_list = []
                pattern = [1]*length
                name = '../%s/%s_%d.pdb'%(dir_name,out_filename, molecule_number)
                for position, character in enumerate(string):
                    if character == 'S':
                        if position == 0:
                            seq_list.append(np.random.choice([stf,strf], p=[0.5, 0.5]))
                        elif position == (length-1):
                            seq_list.append(np.random.choice([stl,strl], p=[0.5, 0.5]))
                        else:
                            seq_list.append(np.random.choice([st,stre], p=[0.5, 0.5]))
                    elif character == 'D':
                        if position == 0:
                            seq_list.append(np.random.choice([dbf,dbrf], p=[0.5, 0.5]))
                        elif position == (length-1):
                            seq_list.append(np.random.choice([dbl,dbrl], p=[0.5, 0.5]))
                        else:
                            seq_list.append(np.random.choice([db,dbre], p=[0.5, 0.5]))
                    elif character == 'M':
                        if position == 0:
                            if protstate == 1:
                                seq_list.append(np.random.choice([ma2f,malf,marf], p=[0.7, 0.15, 0.15]))
                            elif protstate == 2:
                                seq_list.append(np.random.choice([ma2f,malf,marf], p=[0.5, 0.25, 0.25]))
                            elif protstate == 3:
                                seq_list.append(np.random.choice([malf,marf], p=[0.5, 0.5]))
                            elif protstate == 4:
                                seq_list.append(np.random.choice([malf,marf,madf], p=[0.4, 0.4, 0.2]))
                            elif protstate == 5:
                                seq_list.append(np.random.choice([malf,marf,madf], p=[0.15, 0.15, 0.7]))
                            elif protstate == 6:
                                seq_list.append(np.random.choice([malf,marf,madf], p=[0.05, 0.05, 0.9]))
                            elif protstate == 7:
                                seq_list.append(madf)
                        elif position == (length-1):
                            if protstate == 1:
                                seq_list.append(np.random.choice([ma2l,mall,marl], p=[0.7, 0.15, 0.15]))
                            elif protstate == 2:
                                seq_list.append(np.random.choice([ma2l,mall,marl], p=[0.5, 0.25, 0.25]))
                            elif protstate == 3:
                                seq_list.append(np.random.choice([mall,marl], p=[0.5, 0.5]))
                            elif protstate == 4:
                                seq_list.append(np.random.choice([mall,marl,madl], p=[0.4, 0.4, 0.2]))
                            elif protstate == 5:
                                seq_list.append(np.random.choice([mall,marl,madl], p=[0.15, 0.15, 0.7]))
                            elif protstate == 6:
                                seq_list.append(np.random.choice([mall,marl,madl], p=[0.05, 0.05, 0.9]))
                            elif protstate == 7:
                                seq_list.append(madl)
                        else:
                            if protstate == 1:
                                seq_list.append(np.random.choice([ma2,mal,mar], p=[0.7, 0.15, 0.15]))
                            elif protstate == 2:
                                seq_list.append(np.random.choice([ma2,mal,mar], p=[0.5, 0.25, 0.25]))
                            elif protstate == 3:
                                seq_list.append(np.random.choice([mal,mar], p=[0.5, 0.5]))
                            elif protstate == 4:
                                seq_list.append(np.random.choice([mal,mar,mad], p=[0.4, 0.4, 0.2]))
                            elif protstate == 5:
                                seq_list.append(np.random.choice([mal,mar,mad], p=[0.15, 0.15, 0.7]))
                            elif protstate == 6:
                                seq_list.append(np.random.choice([mal,mar,mad], p=[0.05, 0.05, 0.9]))
                            elif protstate == 7:
                                seq_list.append(mad)
                # run the copolymer random walk method
                polymer = copolymer(seq_list, length, pattern=pattern, forcefield=f, density=0.1, limit=0.2, capped=False, unwrap=True)
                # write a pdb file    
                polymer.write_pdb(name)
                print('\n')
                print(name+' generated')
    elif seq_file == 'not_specified':
        if (protstate > 7) or (protstate < 1):
            while True:
                try:
                    print "\n"
                    print "Protonation states of maleic acid residues"
                    print "number    1     2     3     4     5     6     7"
                    print "pH        5     6     7     8     9    10   >10"
                    print "charge -0.3  -0.5  -1.0  -1.2  -1.7  -1.9  -2.0"
                    protstate = int(raw_input("Select the number of protonation state: "))
                except ValueError:
                    print('\n')
                    print("Sorry, I didn't understand that.")
                    continue
                if (protstate > 7) or (protstate < 1):
                    print('\n')
                    print("Wrong number. Select the one from the table.")
                    continue
                else:
                    break
        if random == 'no':
            ##############
            ## MANUALLY ##
            ##############
            # You can set sequence here or you'll be asked to type it in terminal:
            
            n = 0
            while True:
                if (n == 0) or not (isinstance(n,int)):
                    try:
                        n = int(raw_input('How many polymer molecules do you want to generate?? Integer number, please. Example: 322\n'))
                    except:
                        continue
                    continue
                else:
                    break
                    
            mol_list = []
            for molecule_number in range(n):
                sequence = 'thisstringdoesntmeananything'
                while True:
                    if len(set(sequence).difference("MSD")) != 0:
                        sequence = raw_input('Enter your sequence. Only "M","S","D" monomers are allowed. Example: MSMMSDMS\n')
                        continue
                    else:
                        break
                mol_list.append(sequence)
                
            for molecule_number, sequence in enumerate(mol_list):
                length = len(sequence)
                seq_list = []
                pattern = [1]*length
                name = '../%s/%s_%d.pdb'%(dir_name,out_filename, molecule_number)        
                for position, character in enumerate(sequence):
                    if character == 'S':
                        if position == 0:
                            seq_list.append(np.random.choice([stf,strf], p=[0.5, 0.5]))
                        elif position == (length-1):
                            seq_list.append(np.random.choice([stl,strl], p=[0.5, 0.5]))
                        else:
                            seq_list.append(np.random.choice([st,stre], p=[0.5, 0.5]))
                    elif character == 'D':
                        if position == 0:
                            seq_list.append(np.random.choice([dbf,dbrf], p=[0.5, 0.5]))
                        elif position == (length-1):
                            seq_list.append(np.random.choice([dbl,dbrl], p=[0.5, 0.5]))
                        else:
                            seq_list.append(np.random.choice([db,dbre], p=[0.5, 0.5]))
                    elif character == 'M':
                        if position == 0:
                            if protstate == 1:
                                seq_list.append(np.random.choice([ma2f,malf,marf], p=[0.7, 0.15, 0.15]))
                            elif protstate == 2:
                                seq_list.append(np.random.choice([ma2f,malf,marf], p=[0.5, 0.25, 0.25]))
                            elif protstate == 3:
                                seq_list.append(np.random.choice([malf,marf], p=[0.5, 0.5]))
                            elif protstate == 4:
                                seq_list.append(np.random.choice([malf,marf,madf], p=[0.4, 0.4, 0.2]))
                            elif protstate == 5:
                                seq_list.append(np.random.choice([malf,marf,madf], p=[0.15, 0.15, 0.7]))
                            elif protstate == 6:
                                seq_list.append(np.random.choice([malf,marf,madf], p=[0.05, 0.05, 0.9]))
                            elif protstate == 7:
                                seq_list.append(madf)
                        elif position == (length-1):
                            if protstate == 1:
                                seq_list.append(np.random.choice([ma2l,mall,marl], p=[0.7, 0.15, 0.15]))
                            elif protstate == 2:
                                seq_list.append(np.random.choice([ma2l,mall,marl], p=[0.5, 0.25, 0.25]))
                            elif protstate == 3:
                                seq_list.append(np.random.choice([mall,marl], p=[0.5, 0.5]))
                            elif protstate == 4:
                                seq_list.append(np.random.choice([mall,marl,madl], p=[0.4, 0.4, 0.2]))
                            elif protstate == 5:
                                seq_list.append(np.random.choice([mall,marl,madl], p=[0.15, 0.15, 0.7]))
                            elif protstate == 6:
                                seq_list.append(np.random.choice([mall,marl,madl], p=[0.05, 0.05, 0.9]))
                            elif protstate == 7:
                                seq_list.append(madl)
                        else:
                            if protstate == 1:
                                seq_list.append(np.random.choice([ma2,mal,mar], p=[0.7, 0.15, 0.15]))
                            elif protstate == 2:
                                seq_list.append(np.random.choice([ma2,mal,mar], p=[0.5, 0.25, 0.25]))
                            elif protstate == 3:
                                seq_list.append(np.random.choice([mal,mar], p=[0.5, 0.5]))
                            elif protstate == 4:
                                seq_list.append(np.random.choice([mal,mar,mad], p=[0.4, 0.4, 0.2]))
                            elif protstate == 5:
                                seq_list.append(np.random.choice([mal,mar,mad], p=[0.15, 0.15, 0.7]))
                            elif protstate == 6:
                                seq_list.append(np.random.choice([mal,mar,mad], p=[0.05, 0.05, 0.9]))
                            elif protstate == 7:
                                seq_list.append(mad)

                # run the copolymer random walk method
                polymer = copolymer(seq_list, length, pattern=pattern, forcefield=f, density=0.1, limit=0.2, capped=False, unwrap=True)
                # write a pdb file    
                polymer.write_pdb(name)
                print('\n')
                print(name+' generated')
        elif random=='yes':
            ############
            ## RANDOM ##
            ############
            # the length of polymer, seq_list with sequence of monomers and pattern which is always [1, 1, 1, ...] are generated
            length = int(np.random.normal(meanlen,rmsd))
            print('\n')
            print('Length of polymer molecule is '+str(length)+'\n')
            seq_list = []
            pattern = [1]*length
            
            if protstate == 1:
                seq_list.append(np.random.choice([ma2f,malf,marf,stf,strf,dbf,dbrf], p=[pm*0.7, pm*0.15, pm*0.15, ps, ps, pd, pd]))
                for i in xrange(length-2):
                    seq_list.append(np.random.choice([ma2,mal,mar,st,stre,db,dbre], p=[pm*0.7, pm*0.15, pm*0.15, ps, ps, pd, pd]))
                seq_list.append(np.random.choice([ma2l,mall,marl,stl,strl,dbl,dbrl], p=[pm*0.7, pm*0.15, pm*0.15, ps, ps, pd, pd]))
            elif protstate == 2:
                seq_list.append(np.random.choice([ma2f,malf,marf,stf,strf,dbf,dbrf], p=[pm*0.5, pm*0.25, pm*0.25, ps, ps, pd, pd]))
                for i in xrange(length-2):
                    seq_list.append(np.random.choice([ma2,mal,mar,st,stre,db,dbre], p=[pm*0.5, pm*0.25, pm*0.25, ps, ps, pd, pd]))
                seq_list.append(np.random.choice([ma2l,mall,marl,stl,strl,dbl,dbrl], p=[pm*0.5, pm*0.25, pm*0.25, ps, ps, pd, pd]))
            elif protstate == 3:
                seq_list.append(np.random.choice([malf,marf,stf,strf,dbf,dbrf], p=[pm*0.5, pm*0.5, ps, ps, pd, pd]))
                for i in xrange(length-2):
                    seq_list.append(np.random.choice([mal,mar,st,stre,db,dbre], p=[pm*0.5, pm*0.5, ps, ps, pd, pd]))
                seq_list.append(np.random.choice([mall,marl,stl,strl,dbl,dbrl], p=[pm*0.5, pm*0.5, ps, ps, pd, pd]))
            elif protstate == 4:
                seq_list.append(np.random.choice([malf,marf,madf,stf,strf,dbf,dbrf], p=[pm*0.4, pm*0.4, pm*0.2, ps, ps, pd, pd]))
                for i in xrange(length-2):
                    seq_list.append(np.random.choice([mal,mar,mad,st,stre,db,dbre], p=[pm*0.4, pm*0.4, pm*0.2, ps, ps, pd, pd]))
                seq_list.append(np.random.choice([mall,marl,madl,stl,strl,dbl,dbrl], p=[pm*0.4, pm*0.4, pm*0.2, ps, ps, pd, pd]))
            elif protstate == 5:
                seq_list.append(np.random.choice([malf,marf,madf,stf,strf,dbf,dbrf], p=[pm*0.15, pm*0.15, pm*0.7, ps, ps, pd, pd]))
                for i in xrange(length-2):
                    seq_list.append(np.random.choice([mal,mar,mad,st,stre,db,dbre], p=[pm*0.15, pm*0.15, pm*0.7, ps, ps, pd, pd]))
                seq_list.append(np.random.choice([mall,marl,madl,stl,strl,dbl,dbrl], p=[pm*0.15, pm*0.15, pm*0.7, ps, ps, pd, pd]))
            elif protstate == 6:
                seq_list.append(np.random.choice([malf,marf,madf,stf,strf,dbf,dbrf], p=[pm*0.05, pm*0.05, pm*0.9, ps, ps, pd, pd]))
                for i in xrange(length-2):
                    seq_list.append(np.random.choice([mal,mar,mad,st,stre,db,dbre], p=[pm*0.05, pm*0.05, pm*0.9, ps, ps, pd, pd]))
                seq_list.append(np.random.choice([mall,marl,madl,stl,strl,dbl,dbrl], p=[pm*0.05, pm*0.05, pm*0.9, ps, ps, pd, pd]))
            elif protstate == 7:
                seq_list.append(np.random.choice([madf,stf,strf,dbf,dbrf], p=[pm, ps, ps, pd, pd]))
                for i in xrange(length-2):
                    seq_list.append(np.random.choice([mad,st,stre,db,dbre], p=[pm, ps, ps, pd, pd]))
                seq_list.append(np.random.choice([madl,stl,strl,dbl,dbrl], p=[pm, ps, ps, pd, pd]))
            
            for molecule_number in xrange(n):
                print('\n')
                print('Building polymer molecule %d'%molecule_number)
                name = '../%s/%s_%d.pdb'%(dir_name,out_filename,molecule_number)
                # run the copolymer random walk method
                polymer = copolymer(seq_list, length, pattern=pattern, forcefield=f, density=0.1, limit=0.2, capped=False, unwrap=True)
                # write a pdb file    
                polymer.write_pdb(name)
                print('\n')
                print(name+' generated')

    os.chdir('../')
    shutil.rmtree('tmp')
