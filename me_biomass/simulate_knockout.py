# ==========================================================

# ============================================================
try:
    from qminospy.me1 import ME_NLP1
except:
    pass

import cobrame
from cobrame.solve.algorithms import solve_at_growth_rate, binary_search

import pickle
import os
from os.path import abspath, dirname, relpath
import json
import time
import argparse
from me_biomass.load_model import load_me_model

# ------------------------------------------------------------
# Modules to Manipulate Model
# ------------------------------------------------------------
import cobra

# ************************************************************
# Parameters
# ------------------------------------------------------------
# Bisection parameters
MU_PREC = 1e-13
MU_MIN = 0.
MU_MAX = 3.0

use_soplex=False

simstr = ''

aux_to_ko = {'default': [],
             'pydxn': ['PDX5PS1', 'PDX5PS2'], # PDX5PS in iJO, but unlumped for ME
             'thm': ['THZPSN31'],
             'nac': ['ASPO3', 'ASPO4', 'ASPO5', 'ASPO6'],
             'thf': ['DHFR'],  # actually GCALLD, but seems unlikely
             'met__L': ['HSST'],  # from flexneri 2a
             'pnto__R': ['PANTS'],
             'ribflv': ['RBFSb'],
             'skm': ['DHQTi'],
             'trp__L': ['IGPS'],
             'leu__L': ['IPMD'],  # from DH10b
             'btn': ['ALLTN', 'DBTS'],
             'arg__L': ['OCBT'],
             'phe__L': ['PPNDH'],
             'his__L': ['HISTD'],
             'thr__L': ['HSK', 'PTHRpp'],
             'asn__L': ['ASNS2', 'ASNS1'],
             'tyr__L': ['PPND'],
             'gln__L': ['GLNS'],
             'glu__L': ['GLUDy', 'GLUSy']
             }

def str2bool(input):
    if input.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif input.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# ************************************************************
# Argument
parser = argparse.ArgumentParser(description='Simulation parameters.')

parser.add_argument('gene', help='', type=str)
parser.add_argument('aerobicity', help='', type=str)
parser.add_argument('auxotrophy', help='', type=str)
parser.add_argument('media', help='', type=str)
parser.add_argument('source', help='', type=str)

args = parser.parse_args()

# pair must be in form of ko1strain1:ko2strain1-ko1strain2:ko2strain2
GENE_ID = args.gene
AEROBICITY = args.aerobicity
AUXOTROPHY = args.auxotrophy
MEDIA = args.media
SOURCE = args.source

# ************************************************************
# Load and prepare model for simulations
here = dirname(abspath(__file__))


model = load_me_model(json=True)
if GENE_ID.lower() == 'none':
    GENE_ID = None

if GENE_ID:
    model.reactions.get_by_id('translation_' + GENE_ID).knock_out()
    print('loaded model w/ knockout')

# Close default source of nutrient
if SOURCE == 'N':
    source_rxn = model.reactions.EX_nh4_e
elif SOURCE == 'C':
    source_rxn = model.reactions.EX_glc__D_e
elif SOURCE == 'P':
    source_rxn = model.reactions.EX_pi_e
elif SOURCE == 'S':
    source_rxn = model.reactions.EX_so4_e
else:
    raise Exception('bad source')
source_rxn.lower_bound = 0
model.reactions.get_by_id(MEDIA).lower_bound = -1000



# ===================Set auxotrophy=========================================
if AUXOTROPHY == 'thr__L':
    model.reactions.EX_gly_e.upper_bound = 0

for r in aux_to_ko[AUXOTROPHY]:
    for rxn in model.process_data.get_by_id(r).parent_reactions:
        print('knocked out', rxn, 'for', AUXOTROPHY)
        rxn.knock_out()

if AUXOTROPHY != 'default':
    # add exchange for met
    model.add_reaction(cobrame.MEReaction('EX_%s_c' % AUXOTROPHY))
    model.reactions.get_by_id('EX_%s_c' % AUXOTROPHY).add_metabolites(
        {'%s_c' % AUXOTROPHY: -1})
    model.objective = model.reactions.get_by_id('EX_%s_c' % AUXOTROPHY)

    if AUXOTROPHY == 'thf':
        r = cobra.Reaction('SK_dhf_c')
        model.add_reaction(r)
        r.add_metabolites({'dhf_c': -1})

    try:
        model.reactions.get_by_id('EX_%s_e' % AUXOTROPHY).upper_bound = 0
    except:
        print('No extracellular exchange for ', AUXOTROPHY)

    aux_uptake_r = model.reactions.get_by_id('EX_%s_c' % AUXOTROPHY)
    aux_uptake_r.lower_bound = -1000

# =================Aerobicity===========================================

if AEROBICITY == 'aerobic':
    pass
elif AEROBICITY == 'anaerobic':
    model.reactions.EX_o2_e.lower_bound = 0
else:
    raise Exception('Bad aerobicity')

# ============================================================
tic = time.time()


if not use_soplex:
    def solve_model(me_nlp, hs):
        # Re-compile expressions with new keffs in S and solve
        return me_nlp.bisectmu(precision=MU_PREC, mumin=MU_MIN, mumax=MU_MAX,
                               basis=hs)


    with open('test_model.pickle', 'wb') as f:
        pickle.dump(model, f)

    me_nlp = ME_NLP1(model, growth_key='mu')
    me_nlp.compiled_expressions = me_nlp.compile_expressions()

    x, status, hs = me_nlp.solvelp(.01)

    muopt, hs, xopt, cache = solve_model(me_nlp, hs)

else:
    sol = solve_at_growth_rate(model, .05)
    if model.solution.status == 'optimal':
        binary_search(model, min_mu=.05, max_mu=MU_MAX, mu_accuracy=MU_PREC)
    else:
        model.solution = None

output_dir = os.getcwd() + '/media_sims/'
if GENE_ID:
    output_file = GENE_ID
else:
    output_file = AEROBICITY + '_' + SOURCE + '_' + MEDIA


if model.solution is not None:
    with open(output_dir + '/%s_sol.json' % output_file, 'w') as f:
        json.dump(model.solution.x_dict, f)

else:
    with open(output_dir + '/%s_sol.json' % output_file, 'w') as f:
        json.dump({}, f)
