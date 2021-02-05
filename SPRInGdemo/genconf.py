#------------------------------------------------------------------
# SPRInG-Simple Polydisperse Residue Input Generator
# Ver: Dec-07-2020
# Author: Vaidyanathan Sethuraman
# To generate the initial configuration file for polydisperse chains
# Use NAMD and LigninBuilder to run the script
# Auxiliary functions are given in make_genpsf.py

# 'None' is a reserved keyword- DONT USE IT for PDB/PSF filenames.
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
from make_genpsf import * # function definitions
#------------------------------------------------------------------

# Read input file
if len(sys.argv) != 2:
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()
print('*****************SPRInG_V1.0****************************')
print('Input file name: ', sys.argv[1])
#------------------------------------------------------------------

# Set defaults
graft_opt = [0]; 
input_pdb = 'none'; input_namd = 'none'; input_prm = 'none'
input_pres = 'none'; input_pp = 'none'; input_lbd = 'none'
itertype = 'single'
def_res = 'none'; seg_name = 'S'; res_terminator = 'none'
casenum,mono_deg_poly,num_chains,fpdbflag,ftopflag,fresflag,fpatflag,\
fpres_constraint,fpp_constraint,disperflag,makepdifile,fnamdflag,\
pmolflag,cleanslate,flbdflag,packtol,maxatt,conftol = def_vals()
#------------------------------------------------------------------

# Read from file
with open(sys.argv[1]) as farg:
    for line in farg:
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        words = line.split()
        # call all functions
        if words[0] == 'case_num': 
            casenum = int(words[1])
        elif words[0] == 'biomass_type': 
            biomas_typ = words[1]
        elif words[0] == 'disperse':
            disperflag = 1
            if words[1] == 'READ':
                disper_fyle = words[2]; 
            elif words[1] == 'CREATE':
                makepdifile = 1; 
                if len(words) < 5:
                    exit('Not enough arguments for PDI: '+line)
                inp_pdival = float(words[2])
                disper_fyle = words[3]
                npdiatt = int(words[4])
                pditolval = 0
                if len(words) == 6:
                    pditolval = float(words[5])
            else:
                exit('ERR: Unknown PDI option: ' + str(line))
        elif words[0] == 'num_resids':
            mono_deg_poly = int(words[1])
        elif words[0] == 'num_chains':
            num_chains = int(words[1])
        elif words[0] == 'seg_name':
            seg_name = words[1]
        elif words[0] == 'grafting':
            if len(words) < 4 or (len(words)-2)%2 != 0:
                exit('Unknown number of graft options: ' + line)
            else:
                graft_opt[0] = int(words[1])
                for wcnt in range(len(words)-2):
                    graft_opt.append(words[wcnt+2])
        elif words[0] == 'tol':
            conftol = float(words[1])
        elif words[0] == 'nattempts':
            maxatt = int(words[1])
        elif words[0] == 'op_style':
            itertype  = words[1]
            if itertype == 'multi':
                iterinc = int(words[2]) if len(words) == 3 \
                          else exit('Args for multi not found: '+line)
        elif words[0] == 'patch_res_constraint':
            fpres_constraint = 1
            input_pres = words[1] 
        elif words[0] == 'patch_patch_constraint':
            fpp_constraint = 1
            input_pp = words[1] 
        elif words[0] == 'pdb_ipfile':
            if len(words) != 3:
                exit('Unknown number of args: ' + line)
            else:
                input_pdb = words[1]; def_res = words[2]
                fpdbflag = 1
        elif words[0] == 'top_ipfile':
            input_top = words[1]; ftopflag = 1
        elif words[0] == 'resid_inp':
            resinpfyle = words[1]; fresflag = 1
        elif words[0] == 'patch_inp':
            patinpfyle = words[1]; fpatflag = 1
        elif words[0] == 'LigninBuilder':
            input_lbd = words[1]; flbdflag = 1
        elif words[0] == 'namd_inp':
            if len(words) != 3:
                exit('Unknown number of arguments: ' + line)
            else:
                input_namd = words[1]; input_prm = words[2]
                fnamdflag = 1
        elif words[0] == 'clean_directories':
            cleanslate = 1 if words[1] == 'Y' \
                      else 0 if words [1] == 'N' \
                           else exit('Unknown args: '+line)
        elif words[0] == 'gen_packmol':
            pmolflag = 1
            if len(words) != 3 and len(words) != 9:
                exit('Unknown number of arguments: '+ line)
            input_packmol = words[1]; trans_list = []
            if words[1] != 'DEF':
                packtol = float(words[2])
            if len(words) > 3:
                for k in range(6):
                    trans_list.append(words[3+k])            
        elif words[0] == 'terminator':
            res_terminator = words[1]
        else:
            exit('Unknown keyword ' + str(words[0]))
#----------------------------------------------------------------------
# Basic flag checks
outflag = check_all_flags(casenum,fresflag,fpatflag,disperflag,\
                          mono_deg_poly,num_chains,fnamdflag,fpdbflag,\
                          ftopflag)
if outflag == -1:
    exit()
#------------------------------------------------------------------

# Output file names (will be generated automatically)
reslist_fname = 'reslist_' + str(casenum) + '.tcl' #all res list
patch_fname = 'patchlist_' + str(casenum)+ '.tcl' #all patch list
log_fname = 'log_' + str(casenum) + '.txt' #log file
#------------------------------------------------------------------

# Get directory info
srcdir = os.getcwd()
head_outdir = srcdir + str('/casenum_') + str(casenum) # main outdir
#------------------------------------------------------------------

# Create main directories and copy required files
if not os.path.isdir(head_outdir):
    os.mkdir(head_outdir)
elif cleanslate:
    print('Removing old: ', head_outdir)
    srcdir
    shutil.rmtree(head_outdir)
    os.mkdir(head_outdir)
print('Beginning chain generation..')
print('Polymer type/Casenumber: %s, %d' %(biomas_typ,casenum))
#------------------------------------------------------------------

# Check initial and pdb file defaults and copy files
allinitflags = find_init_files(fpres_constraint,fpp_constraint,\
                               fpdbflag,fnamdflag,flbdflag,\
                               makepdifile,input_top,input_pdb,\
                               input_pres,input_pp,input_lbd)
if allinitflags == -1:
    exit()
gencpy(srcdir,head_outdir,sys.argv[1])
gencpy(srcdir,head_outdir,resinpfyle)
gencpy(srcdir,head_outdir,patinpfyle)
gencpy(srcdir,head_outdir,input_top)
if fpres_constraint == 1:
    gencpy(srcdir,head_outdir,input_pres)
if fpp_constraint == 2:
    gencpy(srcdir,head_outdir,input_pp)
if flbdflag == 1:
    gencpy(srcdir,head_outdir,input_lbd)
#------------------------------------------------------------------

# Make monomer array for all chains
if disperflag:
    if makepdifile == 1:
        print("Making chains with target polydispersity...")
        init_pdi_write(inp_pdival,mono_deg_poly,num_chains,disper_fyle\
                       ,npdiatt,pditolval)
        pdigenflag = compile_and_run_pdi(head_outdir)
        if pdigenflag == -1:
            exit()

    print('Polydispersity file: ', disper_fyle)
    deg_poly_all,pdival = make_polydisp_resids(disper_fyle,num_chains)
    if pdival == 0:
        exit()
    gencpy(srcdir,head_outdir,disper_fyle) # copy files to headdir

else:
    deg_poly_all = [mono_deg_poly]*num_chains
    pdival = 1.0
    print('Monodispersed case')
print('Tot ch/res/pat/pdi',num_chains,sum(deg_poly_all),\
      sum(deg_poly_all)-num_chains,pdival)
#------------------------------------------------------------------

# Open log file
flog = open(head_outdir + '/' + log_fname,'w')
init_logwrite(flog,casenum,biomas_typ,deg_poly_all,input_top,\
              seg_name,num_chains,maxatt,conftol,itertype,\
              fpres_constraint,fpp_constraint,resinpfyle,patinpfyle\
              ,disperflag,pdival)
#------------------------------------------------------------------

# Check NAMD inputs and pdb file checks
if fnamdflag == 1:
    if not os.path.exists(input_namd) or not os.path.exists(input_prm):
        exit('No NAMD file found \n')
if fpdbflag and fnamdflag: # if initial pdb file is present
    pdbfyleflag = check_pdb_defaults(input_pdb,def_res,seg_name)
    if pdbfyleflag == -1:
        exit()
    flog.write('NAMD runs will be added in the script\n')
    flog.write('Input PDB name: %s\n' %(input_pdb))
    flog.write('Segment name: %s\n' %(segname))
    gencpy(srcdir,head_outdir,input_pdb) #copy initial pdbfile
#------------------------------------------------------------------

# residue % dictionary mode
resperc_dict = residue_ratios(resinpfyle) 
if not bool(resperc_dict):
    exit('ERROR: Check input residue file \n')
if graft_opt[0] == 1:
    flog.write('Building branched chains...\n')
else:
    flog.write('Building linear chains...\n')
#------------------------------------------------------------------

# patches %dict mode
patchperc_dict = patch_ratios(graft_opt,resperc_dict,patinpfyle) 
if not bool(patchperc_dict):
    exit('ERROR: Check input patch file \n')
#------------------------------------------------------------------

# Open file and write headers
fresin = open(head_outdir + '/' + reslist_fname,'w')
fresin.write(';# Contains all segments for NAMD files.\n')
fpatchin = open(head_outdir + '/' + patch_fname,'w')
fpatchin.write(';# Contains all patches for NAMD files.\n')
#------------------------------------------------------------------
            
# Create cumulative probability distribution of segments/patches
flog.write('Making cumulative distribution for segments..\n')
print('Making cumulative distribution for segments..')
cumul_resarr = cumul_probdist(resperc_dict,flog)
flog.write(str(cumul_resarr)+'\n')

flog.write('Making cumulative distribution for patches..\n')
print('Making cumulative distribution for patches..')
cumul_patcharr = cumul_probdist(patchperc_dict,flog)
flog.write(str(cumul_patcharr)+'\n')
#------------------------------------------------------------------
    
# Set 2D default list and generate segments/patches
res_list = [[] for i in range(num_chains)]
patch_list = [[] for i in range(num_chains-1)]
#------------------------------------------------------------------

# Create residues and check for avg probability 
print('Generating residues..')
flog.write('Creating residue list..\n')
res_list = create_residues(fresin,deg_poly_all,num_chains,seg_name,\
                           resperc_dict,cumul_resarr,conftol,maxatt,\
                           flog,graft_opt,def_res,res_terminator)
if res_list == -1:
    exit()
#------------------------------------------------------------------

# Read patch-patch constraints (in one go)
if fpp_constraint == 1:
    print('Reading patch-patch constraints..')
    flog.write('Reading patch-patch constraints..\n')
    flog.write('All constraints \n')
    ppctr_list = read_patch_incomp(input_pp)
    flog.writelines('\t'.join(str(jval) for jval in ival) +\
                    '\n' for ival in ppctr_list)
else:
    ppctr_list = [[]]
#------------------------------------------------------------------

# Create patches with constraints
print('Generating patches..')
flog.write('Creating patches list..\n')
patch_list = create_patches(fpatchin,deg_poly_all,num_chains,seg_name,\
                            patchperc_dict,cumul_patcharr,conftol,\
                            maxatt,flog,fpres_constraint,fpp_constraint,\
                            input_pres,res_list,ppctr_list,graft_opt)

if patch_list == -1:
    exit()
#------------------------------------------------------------------

# Write to file: Headers
flog.write('Writing data to files \n')
flog.write('Output style %s\n' %(itertype))
print('Writing data to files..')
print('Output style: ', itertype)
#------------------------------------------------------------------

#PACKMOL directives
if pmolflag: 
    flog.write('Initiating packmol files for %d\n..' %(casenum))
    print('Initiating packmol files for ', casenum)
    if input_packmol == 'DEF':
        packname = 'packmol_'+biomas_typ+'_case_'+str(casenum)+'.inp'
    else:
        packname = input_packmol
    fpack = open(head_outdir + '/' + packname,'w')        
    initiate_packmol(fpack,biomas_typ,num_chains,packtol)
#------------------------------------------------------------------

# Make tcl output directory and auxiliary files
tcldir = head_outdir
if not os.path.isdir(tcldir):
    os.mkdir(tcldir)
fbund = make_auxiliary_files(tcldir,biomas_typ,num_chains,input_top,\
                             flbdflag,input_lbd)
if flbdflag == 1:
    gencpy(srcdir,tcldir,input_lbd)
#------------------------------------------------------------------

# Write for each chain
for chcnt in range(num_chains):
    chnum = chcnt + 1
    flog.write('****Writing chain number: %d***\n' %(chnum))
    print('Writing chain number: ', chnum)

    #prefix for pdb/psf/tcl files
    pdbpsf_name = biomas_typ + '_chnum_' + str(chnum) 
    tcl_fname  =  pdbpsf_name +'.tcl' 

    # Copy NAMD files
    if fnamdflag:
        gencpy(srcdir,tcldir,input_prm) 
        fr = open(input_namd,'r')
        fw = open('mini.conf','w')
        fid = fr.read().replace("py_inpname",pdbpsf_name)
        fw.write(fid)
        fr.close(); fw.close()
        gencpy(srcdir,tcldir,'mini.conf')

    deg_poly_this_chain = deg_poly_all[chcnt]
    psfgen_headers(fbund,input_top,pdbpsf_name)
    if itertype == 'single':
        fbund.write('%s\t %s\n' %('set outputname', pdbpsf_name))
        
        flog.write('Writing config for  %d chains\n' %(num_chains))
        write_multi_segments(fbund,-1,deg_poly_this_chain,num_chains,chnum,\
                             seg_name,res_list,patch_list,graft_opt,\
                             deg_poly_this_chain)
        psfgen_postprocess(fbund,itertype,0,seg_name,fnamdflag,input_pdb)

    elif itertype == 'multi':
        flog.write('Iteration increment counter %d\n' %(iterinc))
        # Write segments according to iteration number
        iter_num = 1
        nmonsthisiter = iterinc
    
        while nmonsthisiter <= deg_poly_this_chain:
            fbund.write('%s\t %s\n' %('set outputname', pdbpsf_name))
            flog.write('Writing config for n-segments: %d\n' %(nmonsthisiter))
            write_multi_segments(fbund,iter_num,nmonsthisiter,num_chains,\
                                 chnum,seg_name,res_list,patch_list,\
                                 graft_opt,deg_poly_this_chain)
            psfgen_postprocess(fbund,itertype,iter_num,seg_name,\
                               fnamdflag,input_pdb)
            if fnamdflag == 1:
                out_namd = 'mini' + str(iter_num) + '.out'
                run_namd(fbund,'namd2','mini.conf',out_namd)
            iter_num  += 1
            nmonsthisiter = nmonsthisiter + iterinc

        # Write the rest in one go
        if deg_poly_this_chain%iterinc != 0:
            fbund.write('%s\t %s\n' %('set outputname', pdbpsf_name))
            flog.write('Writing config for n-segments: %d\n' \
                       %(deg_poly_this_chain))
            iter_num += 1
            write_multi_segments(fbund,iter_num,deg_poly_this_chain,\
                                 num_chains,chnum,seg_name,res_list,\
                                 patch_list,graft_opt,deg_poly_this_chain)
            psfgen_postprocess(fbund,itertype,iter_num,seg_name,\
                               fnamdflag,input_pdb)

            if fnamdflag == 1:
                out_namd = 'mini' + str(iter_num) + '.out'
                run_namd(fbund,'namd2','mini.conf',out_namd)

    else:
        exit('ERROR: Unknown output write style option: ' + itertype)
  
    if pmolflag:
        make_packmol(fpack,pdbpsf_name,1,trans_list)
#------------------------------------------------------------------
if flbdflag == 1:
    fbund.write(';# Using LigninBuilder to generate init config\n')
    fbund.write('package require ligninbuilder\n')
    fbund.write('::ligninbuilder::makelignincoordinates . . \n')

fbund.write('exit\n')
fbund.close() #Exit and close file
#------------------------------------------------------------------

#Extra PACKMOL directives
if pmolflag:
    fpack.write('\n')
    fpack.write('#---End of PACKMOL file----\n')
    fpack.close()
flog.write('Completed psf generation for casenum: %d\n' %(casenum))
print('Completed psf generation for casenum: ', casenum)
#------------------------------------------------------------------

# Close files
fresin.close()
fpatchin.close()
flog.close()
#------------------------------------------------------------------
