#!/usr/bin/env python
# from __future__ import print_function
import parmed as pmd
from weighted_levenshtein import lev
from parmed.charmm import *
import glob
import numpy as np

# Non default cost matricies for levenshtein distances.
delcost = np.ones(128, dtype=np.float64)
subcost = np.ones((128, 128), dtype=np.float64)
delcost[ord('L')] = 0.1  # Removing an L to find a match should be low cost
ellist = [ord('C'), ord('O'), ord('H')]

for el1 in ellist:
    for el2 in ellist:
        subcost[el1, el2] = 100  # Don't allow substitutions between elements.


def findmissingparameters(self, parmset):
    """
        Based on ParmEd's "load_parameters" routine, this goes through the loaded psf and determines what 
        parameters are missing.
        The list of missing parameters is passed back to the caller
        """
    missingatomprms = set()
    # First load the atom types
    for atom in self.atoms:
        try:
            if isinstance(atom.type, int):
                atype = parmset.atom_types_int[atom.type]
            else:
                atype = parmset.atom_types_str[atom.type.upper()]
        except KeyError:
            missingatomprms.add(atom.type)
    missingbondprms = set()
    # Next load all of the bonds
    for bond in self.bonds:
        # Construct the key
        key = (min(bond.atom1.type, bond.atom2.type),
               max(bond.atom1.type, bond.atom2.type))
        try:
            bond.type = parmset.bond_types[key]
        except KeyError:
            missingbondprms.add(key)
    # Next load all of the angles. If a Urey-Bradley term is defined for
    # this angle, also build the urey_bradley and urey_bradley_type lists
    missingangleprms = set()
    for ang in self.angles:
        # Construct the key
        key = (min(ang.atom1.type, ang.atom3.type), ang.atom2.type,
               max(ang.atom1.type, ang.atom3.type))
        try:
            ang.type = parmset.angle_types[key]
        except KeyError:
            missingangleprms.add(key)
    # Next load all of the dihedrals.
    missingdihedralprms = set()
    for dih in self.dihedrals:
        # Store the atoms
        a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
        key = (a1.type, a2.type, a3.type, a4.type)
        # First see if the exact dihedral is specified
        if (not key in parmset.dihedral_types) and (not key[::-1] in parmset.dihedral_types):
            # Check for wild-cards
            oldkey = key
            key = ('X', a2.type, a3.type, 'X')
            if not key in parmset.dihedral_types:
                missingdihedralprms.add(oldkey)
    return [missingatomprms, missingbondprms, missingangleprms, missingdihedralprms]


CharmmPsfFile.findmissingparameters = findmissingparameters


def removeduplicates(inputset):
    # ParmEd does not keep track of ordering. This just makes sure that only truly unique parameters are in the set.
    returnset = set()
    for el in inputset:
        if not (el in returnset or el[::-1] in returnset):
            returnset.add(el)
    return returnset


def missingparameterstext(psflist, parameterlist):
    # Pass a directory that has psfs in it and a list of parameter files that contain analogous parameters, and return 
    # a string containing all the parameters that need to be added for simulation.
    params = CharmmParameterSet()
    for p in parameterlist:
        params.read_parameter_file(p)
    returntext = "! This file was written by analogy from the following input parameter files: {}" \
                 "\n\n".format(str(prmlist)[1:-1])
    missingset = [set(), set(), set(), set()]
    for psf in psflist:
        mol = pmd.load_file(psf)
        newmissing = mol.findmissingparameters(params)
        for i in range(4):
            for el in newmissing[i]:
                missingset[i].add(el)
    # Eliminate duplicate dihedrals.
    if len(missingset[0]):
        print(
            "The following atomtypes are missing nonbonded terms from the input parameter files. This is probably an "
            "input problem, and so we are exiting.")
        print(missingset[0])
        exit()
    for i in range(1, 4):
        missingset[i] = removeduplicates(missingset[i])
    # Write bonds section
    returntext += "BONDS\n"
    for missingbondtype in missingset[1]:
        mindistance = 1000
        typekey = "-".join(missingbondtype)
        for k in params.bond_types:
            compkey = "-".join(k)
            ed = lev(typekey, compkey, substitute_costs=subcost, delete_costs=delcost)
            if ed < mindistance:
                replacementtype = k
                mindistance = ed
        returntext += "%-8s %-8s %.3f %.4f ! From %-8s %-8s\n" % (
        missingbondtype[0], missingbondtype[1], params.bond_types[replacementtype].k,
        params.bond_types[replacementtype].req, replacementtype[0], replacementtype[1])
    # Write angles section
    returntext += "\nANGLES\n"
    for missingangletype in missingset[2]:
        mindistance = 1000
        typekey = "-".join(missingangletype)
        for k in params.angle_types:
            compkey = "-".join(k)
            ed = lev(typekey, compkey, substitute_costs=subcost, delete_costs=delcost)
            if ed < mindistance:
                replacementtype = k
                mindistance = ed
        if params.urey_bradley_types[replacementtype].k == 0:
            returntext += "%-8s %-8s %-8s %.3f %.4f ! From %-8s %-8s %-8s\n" % (
            missingangletype[0], missingangletype[1], missingangletype[2], params.angle_types[replacementtype].k,
            params.angle_types[replacementtype].theteq, replacementtype[0], replacementtype[1], replacementtype[2])
        else:
            returntext += "%-8s %-8s %-8s %.3f %.4f %.2f %.4f ! From %-8s %-8s %-8s\n" % (
            missingangletype[0], missingangletype[1], missingangletype[2], params.angle_types[replacementtype].k,
            params.angle_types[replacementtype].theteq, params.urey_bradley_types[replacementtype].k,
            params.urey_bradley_types[replacementtype].req, replacementtype[0], replacementtype[1], replacementtype[2])
    # Write dihedrals section
    returntext += "\nDIHEDRALS\n"
    for missingangletype in missingset[3]:
        mindistance = 1000
        typekey = "-".join(missingangletype)
        for k in params.dihedral_types:
            compkey = "-".join(k)
            compkey2 = "-".join(k[::-1])
            ed = lev(typekey, compkey, substitute_costs=subcost, delete_costs=delcost)
            ed2 = lev(typekey, compkey2, substitute_costs=subcost, delete_costs=delcost)
            if ed < mindistance or ed2 < mindistance:
                replacementtype = k
                mindistance = ed
        for prm in params.dihedral_types[replacementtype]:
            returntext += "%-8s %-8s %-8s %-8s %.4f %d %5.1f ! From %-8s %-8s %-8s %-8s\n" % (
            missingangletype[0], missingangletype[1], missingangletype[2], missingangletype[3], prm.phi_k, prm.per,
            prm.phase, replacementtype[0], replacementtype[1], replacementtype[2], replacementtype[3])
    return returntext + "\n"


prmlist = ["../par_lignin.prm", "../par_all36_cgenff.prm"]

fout = open("extraparameters.prm", "w")
fout.write(missingparameterstext(glob.glob("./*psf"), prmlist))
fout.close()
