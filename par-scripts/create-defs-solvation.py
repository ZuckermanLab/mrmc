#!/usr/bin/env python

import sys
#to install parmed: python setup.py install --home=~
#To use parmed: echo $PYTHONPATH=~/lib/python:$PYTHONPATH
from parmed.amber import *
import networkx as nx

#import parmed
#recursive. 
#create a dictionary of edges (i,j) such that each edge gives the set of other edges with which it shares an atom
	

parm = LoadParm(sys.argv[1])
resname = sys.argv[2] #residue name for definitions file
typestart = int(sys.argv[3])
defsfname = sys.argv[4] 
classstart = int(sys.argv[5])
solvfname = sys.argv[6]
#dihcutoff = float(sys.argv[5])
natom = parm.ptr('natom')
ntypes = parm.ptr('ntypes')
#need to get all the bonds, angles, and dihedrals
nbond = parm.ptr('nbonh') + parm.ptr('mbona')

molgraph=nx.Graph()
molgraph.add_nodes_from([0,natom])
for ibond in range(0,nbond):
        iatom = parm.bonds[ibond].atom1.idx
        jatom = parm.bonds[ibond].atom2.idx
	molgraph.add_edge(iatom,jatom)

cycles= nx.cycle_basis(molgraph)
ring_set_atoms=set()
ring_set_bonds=set()
print cycles
for cycle in cycles:
	n=len(cycle)
	for i in range(0,n):
		iatom=cycle[i]
		if (i==n-1):
			jatom=cycle[0]
		else:
			jatom=cycle[i+1]
		ring_set_atoms.add(iatom)
		ring_set_bonds.add((iatom,jatom))


#ring_set_atoms, ring_set_bonds = find_cycles(parm,natom)
print ring_set_atoms
print ring_set_bonds

#sys.exit(0)
#to do the definitions file
defsfile = open(defsfname,'w')
defsfile.write('RESI {0}\n'.format(resname))
for iatom in range(0,natom):
	type = iatom + typestart
	defsfile.write('\tATOM {0} {1:d}\n'.format(parm.atoms[iatom].name,type))
for ibond in range(0,nbond):
        iatom = parm.bonds[ibond].atom1         
        jatom = parm.bonds[ibond].atom2
	#determine if the bond is rotatable -- both atoms must be bonded to at least two others and not be part of a ring
	if ((len(iatom.bond_partners)>1) and (len(jatom.bond_partners)>1) and ((iatom.idx,jatom.idx) not in ring_set_bonds) and ((jatom.idx,iatom.idx) not in ring_set_bonds)):
		rotatable = True 
		#changed my mind: every bond not part of the ring should be rotatable (including double bonds)
		#that way we can deal with twisting motions, letting the force field restrain the motion
		#Examine all the dihedrals whose central bond is (iatom,jatom).  If any have k>cutoff, 
		#for dih in iatom.dihedrals:
		#	if (((dih.atom2 is iatom) and (dih.atom3 is jatom)) or ((dih.atom2 is jatom) and (dih.atom3 is iatom))):
		#		if (dih.type.phi_k>dihcutoff):
		#			rotatable = False
		#			break	
	else:
		rotatable = False
	#todo: logic to determine which bonds are rotatable
	if (rotatable):
		defsfile.write('\tBOND {0} {1} ROTATABLE SIDECHAIN\n'.format(iatom.name,jatom.name))
	else:
		defsfile.write('\tBOND {0} {1}\n'.format(iatom.name,jatom.name))
defsfile.write('END\n')
defsfile.close()


#create solvation file containing EEF1-style hydration volumes for Garden and Zhorov method
#we need to assign class numbers in the same way as convert-to-tinker.py, so this is a duplicat eof the code
#also store the first atom for each class, so as to be able to analyze bonding
nextclass = classstart
classname2num = {}
firstatom = {}
for iatom in range(0,natom):
        #tinkerfile.write(parm.atoms[iatom].name, parm.atoms[iatom].type, parm.atoms[iatom].mass, parm.atoms[iatom].charge, parm.atoms[iatom].nb_idx
        type = iatom + typestart
        classname = parm.atoms[iatom].type
        if (classname not in classname2num):
                iclass = nextclass
                classname2num[classname] = nextclass
                #tinkerfile.write('# {0:d} {1}'.format(iclass,classname)
                nextclass = nextclass + 1
		firstatom[classname] = iatom
        else:
                iclass = classname2num[classname]

print classname2num
print firstatom

#hydration volume and shell thickness
hvol = {}
hshell = {}
for classname in classname2num:
        iclass = classname2num[classname]
	iatom = firstatom[classname]
	atom =  parm.atoms[iatom]
	atnum = atom.atomic_number 
	nbond = molgraph.degree(iatom)
	print iatom,parm[iatom].name,classname,atnum,nbond
	#assume it isn't ionizable for now
	hshell[iclass] = 3.5
	if (atnum==6): 
		#carbon, can't be ionizable
		if (nbond==4):
			#sp3 hybrid carbon, parameters from "CT2" in solvpar22.inp -- possibly check whether it is primary/secondary/tertiary?
			hvol[iclass] = 22.4
		elif (nbond==3):
			#sp2 hybrid carbon.  Need to tell if it's a carbonyl carbon or an aromatic-type carbon.
			bonded = molgraph.neighbors(iatom)
			#assume it's aromatic type -- this would also cover double bonds not in rings, not sure if we really need a different volume for this
			carbonyl = False
			hvol[iclass] = 18.4 
			for jatom in bonded:
				if ((parm.atoms[jatom].atomic_number==8) and (molgraph.degree(jatom)==1)):
					carbonyl = True
					hvol[iclass] = 14.7
					break
		else:
			print "does this structure have sp hybrid carbons or some other invalidity?"
			sys.exit(-1)
	elif (atnum==7):
		#nitrogen
		if (nbond==4):
			#it must be ionized!
			hshell[iclass] = 6.0
			hvol[iclass] = 11.2 #from type NH3 in solvpar22.inp
		elif (nbond==3):
			#it could be an sp3 hybrid nitrogen that is not ionized.  But all sp3 hybrid nitrogens have pKa's greater than 7, so are ionized at pH 7.  
			#therefore we assume it's an sp2 hybrid.  We should determine if it's a guanidinium/imino-type sp2 hybrid or aromatic-type (count the number of carbons)
			bonded = molgraph.neighbors(iatom)
			ccount = 0
			for jatom in bonded:
				if (parm.atoms[jatom].atomic_number==6):
					ccount = ccount + 1
			if (ccount>1):
				#it's an aromatic nitrogen, copy NR1 and NR2 in solvpar22.inp
				hvol[iclass] = 4.4
			else:
				#it's a guanidinium-type
				hvol[iclass] = 11.2
				hshell[iclass] = 6.0
		else:
			print "sp hybrid nitrogen?"
			sys.exit(-1)
	elif (atnum==8):
		#oxygen
		hvol[iclass] = 10.8
		#need to determine if this is part of a carboxyl group.  Look to see if we are bonded to a carbon, then count the oxygens bonded to it.
		if (nbond==1):
			jatom=molgraph.neighbors(iatom)[0]
			if (parm.atoms[jatom].atomic_number==6):
				bonded=molgraph.neighbors(jatom)
				ocount=0
				for katom in bonded:
					if (parm.atoms[jatom].atomic_number==8):
						ocount = ocount+1
				if (ocount==2):
					#it's a carboxyl
					hshell[iclass]=6.0
	elif (atnum==9):
		#fluorine -- this is not in the EEF1 set, but we need it for fluorinated drugs.  The only possibility is a C-F bond.
		#we assume a vdw radius for fluorine of 1.75, for carbon of 1.908 A,  and a c-f bond length of 1.35 A
		#if we use equation 3 in the EEF1 paper we get 16.4 A^3 for fluorine, but this is bigger than oxygen
		#I'm going to pick 8 A^3 for a guess (it should be smaller than oxygen)
		hvol[iclass] = 8.0
	elif (atnum==16):
		#sulfur, see S and SM in solvpar22.inp
		hvol[iclass] = 14.7 
	elif (atnum==1):
		pass #don't bother
	else:
		print "unrecognized element"
		sys.exit(-1)

	

solvfile = open(solvfname,'w')

for classname in classname2num:
        iclass = classname2num[classname]      
	iatom = firstatom[classname]
        atom =  parm.atoms[iatom]
        atnum = atom.atomic_number
	if (atnum>1):
		solvfile.write("{0:d} {1} {2:.1f} {3:.1f}\n".format(iclass,classname,hvol[iclass],hshell[iclass]))

solvfile.close()

