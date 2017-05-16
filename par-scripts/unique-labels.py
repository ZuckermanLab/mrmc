#!/usr/bin/env python

import sys
import networkx as nx
import numpy as np
import math
from networkx.algorithms import isomorphism

#PDB made from the SDF file by babel -- connectivity only. [A more elegant approach would be to use pybel to read the SDF directly, then copy the info into a networkx graph object,
#but I don't have access to pybel.]
pdbfromsdf=sys.argv[1]  
drug=sys.argv[2]
newpdb=sys.argv[3]
equivfname=sys.argv[4]
infile=open(pdbfromsdf, 'r')
#the string indices are zero based, compared to the 1-based column numbers in PDB doucumentation.
#but the ranges do not include the final index
sdfcoords={}
sdfelement={}
sdfcharge={}
atomcount={}
sdfgraph=nx.Graph()
names={}
for line in infile:
	#print line[0:6]
	if ((line[0:6]=='ATOM  ') or (line[0:6]=='HETATM')):
		#it's a set of atomic coordinates:
		index=int(line[6:11])
		sdfgraph.add_node(index)
		el=line[12:16].lstrip().rstrip()
		x=float(line[30:38])
		y=float(line[38:46])
		z=float(line[46:54])
		q=line[78:80]
		print el, x, y, z
		sdfelement[index] = el
		sdfcoords[index] = [x, y, z]
		sdfcharge[index] = q
		if (el!='H'):
			if (el in atomcount):
				atomcount[el] = atomcount[el] + 1
			else:
				atomcount[el] = 1
			s='{0}{1:d}'.format(el,atomcount[el])
			names[index]=s[0:3]
	elif (line[0:6]=='CONECT'):
		i1=int(line[6:11])
		#columns 7-11, 12-16, 17-21, 22-26 may be indices of bonded atoms.
		for start in [11,16,21,26]:
			end=start+5
			s2=line[start:end].lstrip().rstrip()
			if (len(s2)>0):
				i2=int(s2)
				sdfgraph.add_edge(i1,i2)

print atomcount
print names
#sys.exit(-1)

#Make a copy of sdfgraph without the hydrogens in order to match to refgraph.
sdfgraph_noh=sdfgraph.copy()
for i in sdfgraph_noh.nodes():
	if (sdfelement[i][0:1]=='H'):
		sdfgraph_noh.remove_node(i)


#print sdfelement
#print sdfcoords
print sdfgraph.number_of_edges()
print sdfgraph_noh.number_of_edges()
print names


outpdb=open(newpdb,'w')
equivfile=open(equivfname,'w')
equivcounter=1
#this will reassign labels for hydrogens to match the atoms to which they are bonded
pdb_has_hydrogens=False
matcher=isomorphism.GraphMatcher(sdfgraph_noh,sdfgraph_noh)
ok=matcher.is_isomorphic()
if (ok):
	isos=matcher.isomorphisms_iter()
	#The "isos" object appears to be a list of dictionaries, each of which gives a possible matchup between atoms in the SDF and atoms in the PDB, based on indices
	#Each is a possible equivalence between atoms in the SDF and atoms in the PDB. (May not respect chirality)
	for m in isos:
		#Match up the atoms in the SDF to the correct names in the PDB and adding names for hydrogen as well
		newnames={}
		nhydro={}
		for i in sdfgraph.nodes():
			if ((not pdb_has_hydrogens) and (sdfelement[i][0:1]=='H')):
				#it's a hydrogen, but the PDB does not have hydrogens.  Find the name of the atom to which it's bonded. Keep track of how many hydrogens we have bonded to this atom.
				#The logic for naming the hydrogens is taken from replace-labels.awk
				bonded=sdfgraph.neighbors(i)[0]
				if (bonded in nhydro):
					nhydro[bonded]=nhydro[bonded]+1
				else:
					nhydro[bonded]=1
				name=names[m[bonded]]
				hid=' ABCDE'[nhydro[bonded]-1]
				#The new name is H1, etc. if the bonded name was C1,  otherwise it is H plus hte name of the atom. Add the contents of hid (space, A, B, etc.)
				if (sdfelement[bonded][0:1]=='C'):
					hname='H'+name[1:len(name)]+hid
				else:
					hname='H'+name+hid
				#trim to a maximum of four characters
				hname=hname[0:4]
				newnames[i]=hname
			else:
				#either the PDB has hydroges, or this is a heavy atom.  Either way the PDB should have a correct name
				newnames[i]=names[m[i]]	
	
		print newnames
		pdbfmt='HETATM{0:5d} {1:4s} {2:3s}  {3:4d}    {4:8.3f}{5:8.3f}{6:8.3f}{7:6.2f}{8:6.2f}          {9:>2s}{10:2s}\n'
		#The first time, we output the new pdb with HETATM records with the original coordinates 
		if (equivcounter==1):
			firstnames=newnames.copy()
			for i in sdfgraph.nodes():
				line=pdbfmt.format(i,firstnames[i],drug,1,sdfcoords[i][0],sdfcoords[i][1],sdfcoords[i][2],1.0,0.0,sdfelement[i],sdfcharge[i])
				outpdb.write(line)
			for i in sdfgraph.nodes():
				bonded=sdfgraph.neighbors(i)
				line='CONECT{0:5d}'.format(i)
				for b in bonded:
					line=line+'{0:5d}'.format(b)
				line=line+'\n'
				outpdb.write(line)
			outpdb.write('END\n')
			firsttime=False
		#also write matching info file
		equivfile.write('EQUIV {0:d}\n'.format(equivcounter))
		for i in sdfgraph_noh.nodes():
			equivfile.write('{0} {1}\n'.format(firstnames[i],newnames[i]))
		equivcounter += 1
				
else:
	print 'Molecules do not match'
	sys.exit(-1)
