#!/usr/bin/env python

import sys
#to install parmed: python setup.py install --home=~
#To use parmed: echo $PYTHONPATH=~/lib/python:$PYTHONPATH
from parmed.amber import *
#import parmed
parm = LoadParm(sys.argv[1])
typestart = int(sys.argv[2])
classstart = int(sys.argv[3])
tinkerfname = sys.argv[4]
tinkerfile = open(tinkerfname,'w')


natom = parm.ptr('natom')
ntypes = parm.ptr('ntypes')
#need to get all the bonds, angles, and dihedrals
nbond = parm.ptr('nbonh') + parm.ptr('mbona')
nangle = parm.ptr('ntheth') + parm.ptr('mtheta')
ndihedral = parm.ptr('nphih') + parm.ptr('mphia')
#tinkerfile.write(natom, ntypes, nbond

nextclass = classstart
atomlinefmt = 'atom {0:d} {1:d} {2} "{3}"  {4:d} {5:.3f} {6:d}\n'
classname2num = {}
for iatom in range(0,natom):
	#tinkerfile.write(parm.atoms[iatom].name, parm.atoms[iatom].type, parm.atoms[iatom].mass, parm.atoms[iatom].charge, parm.atoms[iatom].nb_idx
	type = iatom + typestart
	classname = parm.atoms[iatom].type
	if (classname not in classname2num):
		iclass = nextclass
		classname2num[classname] = nextclass
		#tinkerfile.write('# {0:d} {1}'.format(iclass,classname)
		nextclass = nextclass + 1
	else:
		iclass = classname2num[classname]
	#todo: the "0" should be the atomic number
	tinkerfile.write(atomlinefmt.format(type, iclass, parm.atoms[iatom].type, parm.atoms[iatom].name, parm.atoms[iatom].atomic_number, parm.atoms[iatom].mass, 1))

for classname in classname2num:
	iclass = classname2num[classname]
	tinkerfile.write('# {0} {1:d}\n'.format(classname,iclass))

bondlinefmt = 'bond {0:d} {1:d} {2:.2f} {3:.4f} #{4} {5}\n'
bondset = set() #to make sure we don't duplicate
for ibond in range(0,nbond):
	iatom = parm.bonds[ibond].atom1
	jatom = parm.bonds[ibond].atom2
	iclass = classname2num[iatom.type]
	jclass = classname2num[jatom.type]
	k = parm.bonds[ibond].type.k
	req = parm.bonds[ibond].type.req
	#todo:really should sort these
	if ((iclass,jclass) not in bondset):
		tinkerfile.write(bondlinefmt.format(iclass,jclass,k,req,iatom.name,jatom.name))
		bondset.add((iclass,jclass))

anglelinefmt = 'angle {0:d} {1:d} {2:d} {3:.2f} {4:.2f} #{5} {6} {7}\n'
angleset = set()	
for iangle in range(0,nangle):
        iatom = parm.angles[iangle].atom1
        jatom = parm.angles[iangle].atom2
	katom = parm.angles[iangle].atom3
        iclass = classname2num[iatom.type]
        jclass = classname2num[jatom.type]
	kclass = classname2num[katom.type]
        k = parm.angles[iangle].type.k
        theteq = parm.angles[iangle].type.theteq
        #todo:really should sort these
        if ((iclass,jclass,kclass) not in angleset):
                tinkerfile.write(anglelinefmt.format(iclass,jclass,kclass,k,theteq,iatom.name,jatom.name,katom.name))
        	angleset.add((iclass,jclass,kclass))

#For the dihedrals, we separate identifying the parameters and writing them out, so that we can come up with multiple dihedrals as needed
ndih = {}
impdih = {}
dihset = set()
for idih in range(0,ndihedral):
	iatom = parm.dihedrals[idih].atom1
        jatom = parm.dihedrals[idih].atom2        
        katom = parm.dihedrals[idih].atom3    
	latom = parm.dihedrals[idih].atom4
	iclass = classname2num[iatom.type]      
        jclass = classname2num[jatom.type]
        kclass = classname2num[katom.type]
	lclass = classname2num[latom.type]
	dihset.add((iclass,jclass,kclass,lclass))
	per = parm.dihedrals[idih].type.per 
	classes = (iclass, jclass, kclass, lclass, per)
	phi_k = parm.dihedrals[idih].type.phi_k
	per = int(parm.dihedrals[idih].type.per)
	phase = parm.dihedrals[idih].type.phase
	dihparams = (phi_k, phase)
	if (parm.dihedrals[idih].improper):
		if (classes not in impdih):
			impdih[classes] = dihparams
	else: #normal dihedral
		if (classes not in ndih):
			ndih[classes] = dihparams

#tinkerfile.write(ndih#
#tinkerfile.write(impdih

imptorsfmt='imptors {0:d} {1:d} {2:d} {3:d} {4:.3f} {5:.1f} {6:d}\n'
#assume there is only one improper dihedral for each combination of classes.
for classes in impdih:
	iclass,jclass,kclass,lclass,per = classes
	phi_k,per = impdih[classes]
	per = int(per)
	tinkerfile.write(imptorsfmt.format(iclass, jclass, kclass, lclass, phi_k, phase, per))
	

#these do not have \n so we can continue
torsfmt1='torsion {0:d} {1:d} {2:d} {3:d} '
torsfmt2=' {0:.3f} {1:.1f} {2:d}' #magnitude, phase, multiplicity
for classes in dihset:
	iclass, jclass, kclass, lclass = classes
	#tinkerfile.write(torsfmt1.format(iclass,jclass,kclass,lclass),
	firsttorsion = True
	for per in range(1,7):
		ndihlookup = (iclass, jclass, kclass, lclass, per)
		if (ndihlookup in ndih):
			phi_k, phase = ndih[ndihlookup]
			if (firsttorsion):
				tinkerfile.write(torsfmt1.format(iclass,jclass,kclass,lclass))
				firsttorsion = False
			tinkerfile.write(torsfmt2.format(phi_k, phase, per))
	if (not firsttorsion):
		tinkerfile.write('\n')
		


vdwfmt="vdw {0:d} {1:.4f} {2:.4f}\n"
#tinker's files have charges to only 4 sig digits
chargefmt="charge {0:d} {1:.4f}\n"

for classname in classname2num:
	classnum = classname2num[classname]
	LJtype = parm.LJ_types[classname]-1
	tinkerfile.write(vdwfmt.format(classnum, parm.LJ_radius[LJtype], parm.LJ_depth[LJtype]))

#Force an integer charge.

totalcharge = 0.0
for itype in range(0,natom):
        type = itype + typestart
        totalcharge += parm.atoms[itype].charge

adjust = -(totalcharge - round(totalcharge,0))/natom
print "Adjusting charges by {0:.6f} to get total charge of {1:.1f}".format(adjust,round(totalcharge,0))

totalcharge2 = 0.0
for itype in range(0,natom):
	type = itype + typestart
        charge =  parm.atoms[itype].charge+adjust
        totalcharge2 += charge
        print parm.atoms[itype].name,charge
	tinkerfile.write(chargefmt.format(type, charge))

print totalcharge2
tinkerfile.close()


