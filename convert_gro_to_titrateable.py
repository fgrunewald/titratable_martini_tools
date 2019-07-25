import re
import random
import argparse
import numpy as np
import MDAnalysis as mda
from numpy.linalg import norm


def u_vect(vect):
    return(vect/norm(vect))

def norm_sphere():
    v_sphere = np.random.normal(0.0, 1, (5000,3))
    return(np.array([ u_vect(vect) for vect in v_sphere]))

def count_prot_res(atoms,resname):
    waters = 0
    for atom in atoms:
        if atom.resname == resname:
           waters = waters + 1
    return(waters)

parser = argparse.ArgumentParser(description='Convert regular to titratable beads')
parser.add_argument('-f',dest='input_file',help='.gro/.trr/.xtc')
parser.add_argument('-o',dest='out_file',help='.gro')
parser.add_argument('-sel', dest='sel',help='name of water atom')
parser.add_argument('-bead',dest='bead_type',help='choice of: water, acid, base, acid-deprot, base-prot')
args = parser.parse_args()

u = mda.Universe(args.input_file)

# 1. Count how many beads need to be converted 
#pH_beads = count_prot_res(u.atoms, args.resname)
pH_beads = len(u.select_atoms(args.sel))

# 2. Do some bookeping

particles={"acid":3,"water":2,"base":2,"acid-deprot":2,"base-prot":3}
total_beads = len(u.atoms) + pH_beads * particles[args.bead_type]

# 3. Define the formats

atom_format='{:>5d}{:<5s}{:>5s}{:5d}{:8.3F}{:8.3F}{:8.3F}{}'


# 4. Write all beads but the protons
selection = u.select_atoms(args.sel)

with open(args.out_file,'w') as out_file:

     protons = []
     atom_count    = 0                
     vectors = norm_sphere()

     out_file.write('titratable MARTINI \n')
     out_file.write('{:>12d}{}'.format(total_beads,'\n'))                      
 
     for atom in u.atoms:
    
         if atom in selection:

            vector = vectors[random.randint(0, len(vectors) - 1)]                                          
            old_pos = 0.1*atom.position                                                                  
            prot = old_pos + vector * 0.17
            dum  = old_pos - vector * 0.17        
            atom_count+=1
            res=atom.resname                         
            out_file.write(atom_format.format(atom.resindex,res,atom.name,atom_count,old_pos[0],old_pos[1],old_pos[2],'\n'))
            atom_count+=1
            out_file.write(atom_format.format(atom.resindex,res,'DN',atom_count,dum[0],dum[1],dum[2],'\n'))
               
            if args.bead_type == 'water':
               protons.append(prot)

            elif args.bead_type == 'acid':
                 v = vector * 0.17
                 v = np.array([v[1],v[0],v[2]])
                 prot_new = old_pos + v
                 protons.append(prot_new)
                 atom_count+=1
                 res=atom.resname   
                 out_file.write(atom_format.format(atom.resindex,res,'DP',atom_count,prot[0],prot[1],prot[2],'\n'))     

            elif args.bead_type == 'acid-deprot':         
                 atom_count+=1
                 res=atom.resname   
                 out_file.write(atom_format.format(atom.resindex,res,'DP',atom_count,prot[0],prot[1],prot[2],'\n'))     

            elif args.bead_type == "base" :
                  atom_count+=1
                  res=atom.resname   
                  out_file.write(atom_format.format(atom.resindex,res,'DP',atom_count,prot[0],prot[1],prot[2],'\n'))     
            
            elif args.bead_type == "base-prot":
                  atom_count += 1
                  res=atom.resname   
                  out_file.write(atom_format.format(atom.resindex,res,'DP',atom_count,prot[0],prot[1],prot[2],'\n'))
                  exc_prot = np.array([-prot[1],prot[0],prot[2]])
                  protons.append(prot)    
         else: 

            pos = 0.1*atom.position
            atom_count+=1
            res=atom.resname   
            out_file.write(atom_format.format(atom.resindex,res,atom.name,atom_count,pos[0],pos[1],pos[2],'\n'))                                                                                                                  
# 5. Write all protons to file         
                        
     resid_count = atom.resindex
     print(atom_count)
     for pos in protons:
#         print('go_here')
         resid_count = resid_count + 1
         atom_count+=1
         out_file.write(atom_format.format(resid_count,'H+','POS',atom_count,pos[0],pos[1],pos[2],'\n'))

# 6. Write box coodinates
     out_file.write('    {:>8.3F}{:>8.3F}{:>8.3F}{}'.format(0.1*u.trajectory[-1].dimensions[0],0.1*u.trajectory[-1].dimensions[1],0.1*u.trajectory[-1].dimensions[2],'\n'))

