#!/usr/bin/env python
from __future__ import print_function

from ligninkmc.kmc_functions import (run_kmc, generate_mol, gen_tcl)
from ligninkmc.create_lignin import (calc_rates, create_initial_monomers, create_initial_events,
                                     create_initial_state, analyze_adj_matrix, adj_analysis_to_stdout)
from ligninkmc.kmc_common import (DEF_E_BARRIER_KCAL_MOL, ADJ_MATRIX, MONO_LIST, MONOMER, OX, GROW, C, Monomer, Event)
import numpy as np

temp = 298.15  # K
rxn_rates = calc_rates(temp, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)

# Set the percentage of S
sg_ratio = 1
pct_s = sg_ratio / (1 + sg_ratio)

# Set the initial and maximum number of monomers to be modeled.
ini_num_monos = 2
max_num_monos = 10

# Maximum time to simulate, in seconds
t_max = 1  # seconds
mono_add_rate = 1e4  # monomers/second

# Use a random number and the given sg_ratio to determine the monolignol types to be initially modeled
monomer_draw = np.random.rand(ini_num_monos)
initial_monomers = create_initial_monomers(pct_s, monomer_draw)

# Initially allow only oxidation events. After they are used to determine the initial state, add 
#     GROW to the events, which allows additional monomers to be added to the reaction at the 
#     specified rate and with the specified ratio
initial_events = create_initial_events(initial_monomers, rxn_rates)
initial_state = create_initial_state(initial_events, initial_monomers)
initial_events.append(Event(GROW, [], rate=mono_add_rate))

result = run_kmc(rxn_rates, initial_state,initial_events, n_max=max_num_monos, t_max=t_max, sg_ratio=sg_ratio)

# Convert the sparse matrix to a full array before printing
print("The adjacency matrix for the simulated lignin is:")
print(result[ADJ_MATRIX].toarray())

# From the list of monomers and the adjacency matrix, we can use LigninKMC to write out a tcl script for psfgen to
# turn into a .psf file.
# fname and sgnames are things that we'd want to change; file name always the same as the segname
gen_tcl(result[ADJ_MATRIX], result[MONO_LIST], toppar_dir="../smilesdemo/toppar/", tcl_fname="psfgen.tcl", psf_fname="L", chain_id="L")  
# Now we switch over into vmd...
