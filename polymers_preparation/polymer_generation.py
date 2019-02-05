from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import copolymer
import numpy as np
import shutil, os, sys
from getopt import getopt

##############################################################################################
#                Script builds styrol-diisobutylene-maleic acid copolymers.
# The @n polymer molecules are constructed with mean length of one molecule @meanlen
# and root mean square deviation of length @rmsd, occurrencies of monomers @ps, @pd and @pm,
# and @protstate from 1 to 7. Output files are .pdb files with name @outfilename_X.pdb.
#  
#   Protonation states of maleic acid residues:
#   number    1     2     3     4     5     6     7
#   pH        5     6     7     8     9    10   >10
#   charge -0.3  -0.5  -1.0  -1.2  -1.7  -1.9  -2.0
#
# @protstate = 0 is an option to create not random polymer. You have to pass @protstate=0 to run() function
# and write sequence of monomers manually. There is an example of how to do it in script.
# @position variable which is passed to monomer build functions declares
# one of three positions monomer can reside: first(@position=1), last(@position=2) and middle(@position=0).
# The difference between positions is that in case of first and last monomers terminal
# hydrogens must not be removed during "polimerization".
# There are 8 types of monomers:
# 4 = (2 different states of protonation) x (2 different carboxyls of maleic acid),
# 2 for diisobutylene monomers with radical group oriented closer to beginning or closer to end of polymer,
# 2 for styrene monomers with benzene group oriented closer to beginning or closer to end of polymer.
##############################################################################################


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

def run(monomers, meanlen, pm, pd, ps, rmsd, name, protstate, i):
    ################################################
    # This function builds one molecule of polymer #
    ################################################

    # NOTE that from this point 'ps' and 'pd' are not the probabilities to have styrol or diisobutylene
    # in any position, but probabilities to have left-oriented or right-oriented styrol\diisobutylene in that position!!
    ps = ps/2
    pd = pd/2



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
    
    # the length of polymer, list with sequence of monomers and pattern which is always [1, 1, 1, ...] are generated
    length = int(np.random.normal(meanlen,rmsd))
    print('\n')
    print('Length of polymer molecule is '+str(length)+'\n')
    list = []
    pattern = [1]*length
    if protstate == 0:
        # This is not a protstate. This is the option to create not random polymers.
        # For example, here 21-monomer-S-and-M-alterating-polymer with 10 M and 11 S (SMSMSMSMSMSMSMSMSMSMS) is contructed.
        # Maleic acids are with charge -1.
        length = 21
        pattern = [1]*length
        list.append(np.random.choice([stf,strf], p=[0.5, 0.5]))
        for i in xrange(length-2):
            if (i%2 == 0):
                list.append(np.random.choice([mal,mar], p=[0.5, 0.5]))
            else:
                list.append(np.random.choice([st, stre], p=[0.5, 0.5]))
        list.append(np.random.choice([stl,strl], p=[0.5, 0.5]))
        ###############
    elif protstate == 1:
        list.append(np.random.choice([ma2f,malf,marf,stf,strf,dbf,dbrf], p=[pm*0.7, pm*0.15, pm*0.15, ps, ps, pd, pd]))
        for i in xrange(length-2):
            list.append(np.random.choice([ma2,mal,mar,st,stre,db,dbre], p=[pm*0.7, pm*0.15, pm*0.15, ps, ps, pd, pd]))
        list.append(np.random.choice([ma2l,mall,marl,stl,strl,dbl,dbrl], p=[pm*0.7, pm*0.15, pm*0.15, ps, ps, pd, pd]))
    elif protstate == 2:
        list.append(np.random.choice([ma2f,malf,marf,stf,strf,dbf,dbrf], p=[pm*0.5, pm*0.25, pm*0.25, ps, ps, pd, pd]))
        for i in xrange(length-2):
            list.append(np.random.choice([ma2,mal,mar,st,stre,db,dbre], p=[pm*0.5, pm*0.25, pm*0.25, ps, ps, pd, pd]))
        list.append(np.random.choice([ma2l,mall,marl,stl,strl,dbl,dbrl], p=[pm*0.5, pm*0.25, pm*0.25, ps, ps, pd, pd]))
    elif protstate == 3:
        list.append(np.random.choice([malf,marf,stf,strf,dbf,dbrf], p=[pm*0.5, pm*0.5, ps, ps, pd, pd]))
        for i in xrange(length-2):
            list.append(np.random.choice([mal,mar,st,stre,db,dbre], p=[pm*0.5, pm*0.5, ps, ps, pd, pd]))
        list.append(np.random.choice([mall,marl,stl,strl,dbl,dbrl], p=[pm*0.5, pm*0.5, ps, ps, pd, pd]))
    elif protstate == 4:
        list.append(np.random.choice([malf,marf,madf,stf,strf,dbf,dbrf], p=[pm*0.4, pm*0.4, pm*0.2, ps, ps, pd, pd]))
        for i in xrange(length-2):
            list.append(np.random.choice([mal,mar,mad,st,stre,db,dbre], p=[pm*0.4, pm*0.4, pm*0.2, ps, ps, pd, pd]))
        list.append(np.random.choice([mall,marl,madl,stl,strl,dbl,dbrl], p=[pm*0.4, pm*0.4, pm*0.2, ps, ps, pd, pd]))
    elif protstate == 5:
        list.append(np.random.choice([malf,marf,madf,stf,strf,dbf,dbrf], p=[pm*0.15, pm*0.15, pm*0.7, ps, ps, pd, pd]))
        for i in xrange(length-2):
            list.append(np.random.choice([mal,mar,mad,st,stre,db,dbre], p=[pm*0.15, pm*0.15, pm*0.7, ps, ps, pd, pd]))
        list.append(np.random.choice([mall,marl,madl,stl,strl,dbl,dbrl], p=[pm*0.15, pm*0.15, pm*0.7, ps, ps, pd, pd]))
    elif protstate == 6:
        list.append(np.random.choice([malf,marf,madf,stf,strf,dbf,dbrf], p=[pm*0.05, pm*0.05, pm*0.9, ps, ps, pd, pd]))
        for i in xrange(length-2):
            list.append(np.random.choice([mal,mar,mad,st,stre,db,dbre], p=[pm*0.05, pm*0.05, pm*0.9, ps, ps, pd, pd]))
        list.append(np.random.choice([mall,marl,madl,stl,strl,dbl,dbrl], p=[pm*0.05, pm*0.05, pm*0.9, ps, ps, pd, pd]))
    elif protstate == 7:
        list.append(np.random.choice([madf,stf,strf,dbf,dbrf], p=[pm, ps, ps, pd, pd]))
        for i in xrange(length-2):
            list.append(np.random.choice([mad,st,stre,db,dbre], p=[pm, ps, ps, pd, pd]))
        list.append(np.random.choice([madl,stl,strl,dbl,dbrl], p=[pm, ps, ps, pd, pd]))
        
    # run the copolymer random walk method
    polymer = copolymer(list, length, pattern=pattern, forcefield=f, density=0.1, limit=0.2, capped=False, unwrap=True)
    
    # write a pdb file    
    polymer.write_pdb(name)
    print('\n')
    print(name+' generated')

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
            
    # During polymerization great amount of data is generated. It makes the main folder
    # look like a Trash. So we build polymer inside tmp folder and then remove it.
    if os.path.isdir('tmp'):
        shutil.rmtree('tmp')
        os.mkdir('tmp')
        os.chdir('tmp')
    else:
        os.mkdir('tmp')
        os.chdir('tmp')
        
    opts, args = getopt(sys.argv[1:], 'n:l:r:p:o:',longopts=['pm=','ps=','pd='])
    # Number of polymer molecules to be generated, their mean length, rmsd of length, protonation state of MA monomers
    n = 100
    meanlen = 36.0
    rmsd = 3.0
    protstate = 3

    # These are the occurrencies of monomers (pm - maleic acid, pd - diisobutylene, ps - styrol)
    pm = 0.25
    pd = 0.75
    ps = 0.0
    out_filename = 'polymer'

    for o, a in opts:
        if o == '-n':
            n = a
        if o == '-l':
            meanlen = float(a)
        if o == '-r':
            rmsd = float(a)
        if o == '-p':
            protstate = a
        if o == '--pm':
            pm = float(a)
        if o == '--ps':
            ps = float(a)
        if o == '--pd':
            pd = float(a)
        if o == '-o':
            out_filename = str(a)
            
       

    print "Normalizing occurrencies..."
    s = pm + ps + pd
    pm /= s
    ps /= s
    pd /= s
    print "Normalized occurrencies are:"
    print('ps = '+str(ps)+'; pm = '+str(pm)+'; pd = '+str(pd))
 

    if (protstate > 7) or (protstate < 0):
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
            
    for i in xrange(int(n)):
        print('\n')
        print('Building polymer molecule %d'%i)
        name = '../polymer_molecules/%s_%d.pdb'%(out_filename,i)
        run(monomers, meanlen=meanlen, pm=pm, pd=pd, ps=ps, rmsd=rmsd, name=name, protstate=protstate, i=i)
    os.chdir('../')
    shutil.rmtree('tmp')
