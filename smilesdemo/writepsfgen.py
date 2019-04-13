from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
monomers = [x for x in Chem.SmilesMolSupplier("monomerdimersmiles/monomers.smiles")]
dimers = [x for x in Chem.SmilesMolSupplier("monomerdimersmiles/dimers.smiles")]
endcaps = [x for x in Chem.SmilesMolSupplier("monomerdimersmiles/endcaps.smiles")]
##This is handy for debugging, mapping indices to the input smiles strings.
# print len(dimers)
# def mol_with_atom_index( mol ):
#     atoms = mol.GetNumAtoms()
#     for idx in range( atoms ):
#         mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx() ) )
#     return mol
# for i in range(len(dimers)):
# 	mol_with_atom_index(dimers[i])
# 	Draw.MolToFile(dimers[i], "test-%d.png" % i)

# for i in range(len(monomers)):
# 	mol_with_atom_index(monomers[i])
# 	Draw.MolToFile(monomers[i], "testm-%d.png" % i)
# for i in range(len(endcaps)):
# 	mol_with_atom_index(endcaps[i])
# 	Draw.MolToFile(endcaps[i], "testend-%d.png" % i)

testmol = [x for x in Chem.SmilesMolSupplier("ligninmol.smiles")]

fout = open("psfgen.tcl", "w")
fout.write("package require psfgen\ntopology toppar/top_all36_cgenff.rtf\ntopology toppar/top_lignin.top\ntopology toppar/top_spirodienone.top\n")
def markC4(mol, name, match):
	if name in ['CAT', 'PHP']:
		atom = mol.GetAtomWithIdx(match[1])
	elif name in ['GUAI']:
		atom = mol.GetAtomWithIdx(match[3])
	elif name in ['SYR', 'TRCN']:
		atom = mol.GetAtomWithIdx(match[9])
	elif name in ['PHPS','GUAS','SYRS']:
		atom = mol.GetAtomWithIdx(match[4])
	if atom.HasProp("residue"):
		atom.SetProp("name", "C4")
def markC1(mol, name, match):
	if name in ['PHP', 'SYR', 'TRCN']:
		atom = mol.GetAtomWithIdx(match[4])
	elif name in ['CAT']:
		atom = mol.GetAtomWithIdx(match[5])
	elif name in ['GUAI', 'GUAS']:
		atom = mol.GetAtomWithIdx(match[7])
	elif name in ['PHPS']:
		atom = mol.GetAtomWithIdx(match[1])
	elif name in ['SYRS']:
		atom = mol.GetAtomWithIdx(match[9])
	if atom.HasProp("residue"):
		atom.SetProp("name", "C1")
#Clearing C1 means that no end-cap modifications are possible.
def clearC1(mol, residue, idxs):
	for idx in idxs:
		atom = mol.GetAtomWithIdx(idx)
		if atom.HasProp("name") and atom.HasProp("residue") and (atom.GetProp("name") == "C1") and (int(atom.GetProp("residue")) == residue):
			#print "Cleared C1 from residue %d" % residue
			atom.ClearProp("name")
for molnum, mol in enumerate(testmol):
	print "Mol: ", molnum
	if len(testmol) > 1:
		segname = "L%d" % molnum
	else:
		segname = "L"
	#Find aromatic units in the test molecule. These are the monomer building blocks we are interested in.
	#Large lignin structures need absurdly high values for maxMatches, since each ring matches with benzene 6 ways.
	benzenerings = mol.GetSubstructMatches(monomers[0], maxMatches=100000)
	#print len(benzenerings)
	for tup in benzenerings:
		for idx in tup:
			atom = mol.GetAtomWithIdx(idx)
			atom.SetUnsignedProp("ringAtom", 1)
	#Find the quinone-type units in the test molecules.
	qrings = mol.GetSubstructMatches(monomers[1], maxMatches=1000)
	for tup in qrings:
		#Skip over the oxygen so that the "ring" still has 6 atoms.
		for idx in tup[1:]:
			atom = mol.GetAtomWithIdx(idx)
			atom.SetUnsignedProp("ringAtom", 1)
	residue = 1
	#Assign residue types and identifiers 
	for monomer in monomers[:0:-1]: #Done in reverse order, since Tricin monomers require special handling.
		matches = mol.GetSubstructMatches(monomer)
		#print matches, monomer.GetProp('_Name')
		if monomer.GetProp('_Name') == 'TRCN':
			#Prune down the "active" part to the part that looks like Syringol
			syringolmatches = mol.GetSubstructMatches(monomers[-2])
			for idx in np.setdiff1d(np.asarray(matches), np.asarray(syringolmatches)):
				atom = mol.GetAtomWithIdx(idx)
				if atom.HasProp("ringAtom"):
					atom.ClearProp("ringAtom")
		#The ordering of the monomers array is very deliberate. If something matches against TRCN, it can't match anything else.
		#Likewise, we compare for SYR, then GUAI, then CAT, then PHP. A monomer can only be assigned once.
		for match in matches:
			newresidue = False
			ringatomcount = 0
			for idx in match:
				atom = mol.GetAtomWithIdx(idx)
				if atom.HasProp("ringAtom"):
					ringatomcount += 1
			#Sometimes, a 4O5 linkage can make G look like S and H look like G. Lets not do that...
			if ringatomcount != 6:
				continue
			for idx in match:
				atom = mol.GetAtomWithIdx(idx)
				if atom.HasProp("ringAtom") and not atom.HasProp("restype"):
					newresidue = True
					atom.SetUnsignedProp("residue", residue)
					atom.SetProp("restype", monomer.GetProp('_Name'))
					atom.SetProp( 'molAtomMapNumber', str( residue ))
			if newresidue:
				markC4(mol, monomer.GetProp('_Name'), match)
				markC1(mol, monomer.GetProp('_Name'), match)
				residue += 1
	#Another debugging draw command.
	#Draw.MolToFile(mol, "testmol-%d.png" % molnum, size=(3600,3600))
	#Write monomers
	fout.write("resetpsf\nsegment %s {\n" % segname)
	residuelist = []
	resnamelist = []
	for atom in mol.GetAtoms():
		if atom.HasProp("ringAtom"):
			residuelist.append(atom.GetProp("residue"))
			resnamelist.append(atom.GetProp("restype"))
	resnamelist = np.array(resnamelist)

	for i, monoletter in enumerate(resnamelist[np.argsort(np.array(residuelist, dtype=np.int))][::6]):
		fout.write("residue %d %s\n" % (i+1, monoletter))
	fout.write("}\n")

	#Now look for dimeric substructures in the lignin.
	for dimer in dimers:
		matches = mol.GetSubstructMatches(dimer)
		#print matches, dimer.GetProp('_Name')
		for match in matches:
			#These are symmetric. Order is irrelevant.
			if dimer.GetProp('_Name') in ['BB', '55']:
				resset = set()
				for idx in match:
					atom = mol.GetAtomWithIdx(idx)
					if atom.HasProp("residue"):
						resset.add(atom.GetProp("residue"))
				resset = list(resset)
				fout.write("patch %s %s:%s %s:%s\n" % (dimer.GetProp('_Name'), segname, resset[0], segname, resset[1]))
				if dimer.GetProp('_Name') == "BB":
					clearC1(mol, int(resset[0]), match)
					clearC1(mol, int(resset[1]), match)
			#For the rest we need to look at the match index. These will tell us what side is residue "1" and which is residue "2" for patching.
			elif dimer.GetProp('_Name') == 'B5':
				atomresid2 = mol.GetAtomWithIdx(match[5]) #The "2" residue
				atomresid1 = mol.GetAtomWithIdx(match[10]) #The "1" residue
				clearC1(mol, int(atomresid1.GetProp("residue")), match)
				oxygenatom = mol.GetAtomWithIdx(match[3])
				oxygenatom.SetUnsignedProp("B5", 1)
				fout.write("patch B5%s %s:%s %s:%s\n" % (atomresid2.GetProp('restype')[0], segname, atomresid1.GetProp('residue'), segname, atomresid2.GetProp('residue')))
			elif dimer.GetProp('_Name') == '4O5':
				atomresida = mol.GetAtomWithIdx(match[5])
				atomresidb = mol.GetAtomWithIdx(match[7])
				if atomresida.HasProp("name") and atomresida.GetProp("name") == "C4" and atomresidb.HasProp("name") and atomresidb.GetProp("name") == "C4":
					print "This should be impossible."
					exit()
				elif atomresida.HasProp("name") and atomresida.GetProp("name") == "C4":
					atomresid1 = atomresida
					atomresid2 = atomresidb
				elif atomresidb.HasProp("name") and atomresidb.GetProp("name") == "C4":
					atomresid2 = atomresida
					atomresid1 = atomresidb
				else:
					print "I should never get here!"
					exit()
				fout.write("patch 4O5 %s:%s %s:%s\n" % (segname, atomresid1.GetProp('residue'), segname, atomresid2.GetProp('residue')))
			elif dimer.GetProp('_Name') in ['AO4', 'AOA', 'GOA']:
				atomresid2 = mol.GetAtomWithIdx(match[6]) #The "2" residue
				atomresid1 = mol.GetAtomWithIdx(match[10]) #The "1" residue
				if dimer.GetProp('_Name') == 'AO4' and mol.GetAtomWithIdx(match[1]).HasProp('B5'):
					continue
				if dimer.GetProp('_Name') == 'GOA':
					fout.write("patch DBAC %s:%s\n" % (segname, atomresid1.GetProp('residue')))
				clearC1(mol, int(atomresid1.GetProp("residue")), match)
				if dimer.GetProp('_Name') != "AO4":
					clearC1(mol, int(atomresid2.GetProp("residue")), match)
				fout.write("patch %s %s:%s %s:%s\n" % (dimer.GetProp('_Name'), segname, atomresid1.GetProp('residue'), segname, atomresid2.GetProp('residue')))
			elif dimer.GetProp('_Name') in ['BO4', 'B1']:
				atomresid1 = mol.GetAtomWithIdx(match[6]) #The "1" residue
				atomresid2 = mol.GetAtomWithIdx(match[10]) #The "2" residue
				clearC1(mol, int(atomresid1.GetProp("residue")), match)
				if dimer.GetProp('_Name') == "B1":
					clearC1(mol, int(atomresid2.GetProp("residue")), match)
				fout.write("patch %s %s:%s %s:%s\n" % (dimer.GetProp('_Name'), segname, atomresid1.GetProp('residue'), segname, atomresid2.GetProp('residue')))
			elif dimer.GetProp('_Name') in ['BDO']:
				atomresid1 = mol.GetAtomWithIdx(match[6]) #The "1" residue
				atomresid2 = mol.GetAtomWithIdx(match[12]) #The "2" residue
				clearC1(mol, int(atomresid1.GetProp("residue")), match)
				fout.write("patch %s %s:%s %s:%s\n" % (dimer.GetProp('_Name'), segname, atomresid1.GetProp('residue'), segname, atomresid2.GetProp('residue')))
			elif dimer.GetProp('_Name') in ['AOG', 'GOG']:
				atomresid1 = mol.GetAtomWithIdx(match[12]) #The "1" residue
				atomresid2 = mol.GetAtomWithIdx(match[6]) #The "2" residue
				if dimer.GetProp('_Name') == 'GOG':
					fout.write("patch DBAC %s:%s\n" % (segname, atomresid1.GetProp('residue')))
				clearC1(mol, int(atomresid1.GetProp("residue")), match)
				clearC1(mol, int(atomresid2.GetProp("residue")), match)
				fout.write("patch %s %s:%s %s:%s\n" % (dimer.GetProp('_Name'), segname, atomresid1.GetProp('residue'), segname, atomresid2.GetProp('residue')))
			elif dimer.GetProp('_Name') in ['SPIR']:
				atomresid1 = mol.GetAtomWithIdx(match[12]) #This is the quinone-type residue
				atomresid2 = mol.GetAtomWithIdx(match[2]) #This is the BO4 linked residue
				atomresid3 = mol.GetAtomWithIdx(match[18]) #This is the funky linked residue.
				clearC1(mol, int(atomresid1.GetProp("residue")), match)
				clearC1(mol, int(atomresid3.GetProp("residue")), match)
				fout.write("patch BO4 %s:%s %s:%s\n" % (segname, atomresid1.GetProp('residue'), segname, atomresid2.GetProp('residue')))
				fout.write("patch SPIR %s:%s %s:%s\n" % (segname, atomresid1.GetProp('residue'), segname, atomresid3.GetProp('residue')))
	#Look for endcaps that aren't the default triol.
	for i, endcap in enumerate(endcaps[::-1]):
		matches = mol.GetSubstructMatches(endcap)
		c1index = 5 - (i/2)
		for match in matches:
			c1atom = mol.GetAtomWithIdx(match[c1index])
			if c1atom.HasProp("name") and c1atom.GetProp("name") == "C1":
				fout.write("patch %s %s:%s\n" % (endcap.GetProp('_Name'), segname, c1atom.GetProp('residue')))
				c1atom.ClearProp("name")
	#All done modifying the psftext through patches. Write it out!
	fout.write("writepsf %s.psf\n" % segname)
fout.close()
