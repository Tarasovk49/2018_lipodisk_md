from sys import argv
from getopt import getopt

opts, args = getopt(argv[1:], 'i:o:')
for o, a in opts:
    if o == '-i':
        in_filename = a
    if o == '-o':
        out_filename = a

def add_ter_to_chains(in_filename,out_filename):
    with open(in_filename) as pdb:
        with open(out_filename,'w') as write:
            resid=1
            for line in pdb.readlines():
                if line[:4]=='ATOM':
                    if int(line[23:26] )<resid:
                        write.write('TER\n')
                        write.write(line)
                        resid = int(line[24:26] )

                    elif int(line[23:26] )==resid:
                        write.write(line)
                    else:
                        resid+=1
                        write.write(line)
                else:
                    write.write(line)

add_ter_to_chains(in_filename,out_filename)
