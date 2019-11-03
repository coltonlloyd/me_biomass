from qminospy.me1 import ME_NLP1
from copy import deepcopy
import pandas as pd
import pickle
import cobra
import numpy as np
from IPython import embed
from me_cofactor.load_model import load_me_model

target_to_flux = {}
target_to_shadow = {}
target_to_reduced = {}
anaerobic = True
if anaerobic:
    suffix = '_anaerobic'
else:
    suffix = ''

model = load_me_model()

if anaerobic:
    model.reactions.EX_o2_e.lower_bound = 0


me_nlp = ME_NLP1(model, growth_key='mu')
me_nlp.compiled_expressions = me_nlp.compile_expressions()
hs = None
for source in ["C", "P", "S", "N"]:
    if source == 'N':
        source_rxn = model.reactions.EX_nh4_e
    elif source == 'C':
        source_rxn = model.reactions.EX_glc__D_e
    elif source == 'P':
        source_rxn = model.reactions.EX_pi_e
    elif source == 'S':
        source_rxn = model.reactions.EX_so4_e
    else:
        raise Exception('bad source')

    source_rxn.lower_bound = 0
    for r in model.reactions.query('EX_'):

        if not r.id.startswith('EX_'):
            continue
        if r.lower_bound < 0:
            print('Skipping and not setting bound of %s to 0' % r.id)
            continue
        met = model.metabolites.get_by_id(r.id.replace('EX_', ''))
        if met.elements.get(source, 0) == 0:
            print('No %s in %s' % (source, met.id))
            continue

        r.lower_bound = -1000
        col = source + suffix + '_' + met.id
        print(col)
        x, status, hs = me_nlp.solvelp(.1, basis=hs)
        if status != 'optimal':
            print(col, 'does not solve', status)
            target_to_flux[col] = {}
            target_to_shadow[col] = {}
            target_to_reduced[col] = {}
            r.lower_bound = 0
            continue

        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-8, mumax=1.5,
                                                 mumin=0, basis=hs)

        target_to_flux[col] = model.solution.x_dict
        target_to_shadow[col] = dict(
            zip(model.metabolites.list_attr('id'), me_nlp.pi))
        target_to_reduced[col] = dict(
            zip(model.reactions.list_attr('id'), me_nlp.rc))
        r.lower_bound = 0

        save_dir = 'substrates%s' % suffix
        pd.DataFrame(target_to_flux).to_csv('./%s/fluxes.csv' % save_dir)
        pd.DataFrame(target_to_shadow).to_csv('./%s/shadow_wrt_target.csv' %
                                              save_dir)
        pd.DataFrame(target_to_reduced).to_csv('./%s/reduced_costs_wrt_target.csv'
                                               % save_dir)
    source_rxn.lower_bound = -1000
