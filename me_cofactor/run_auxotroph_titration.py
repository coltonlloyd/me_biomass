from qminospy.me1 import ME_NLP1
from copy import deepcopy
import pandas as pd
import pickle
import cobra
import numpy as np
from sympy import Basic
from IPython import embed
import cobrame
from multime.auxotroph_analysis.load_model import load_me_model, currency_met_to_synthesis_rxn

me = load_me_model()

aas = [i.replace('_c', '') for i in me.process_data.b2020.amino_acid_count.keys()]
mets = ['quln', 'thf', 'thm', 'fad', 'btn', 'hemed', 'fmn', 'nac', 'thmpp', 'nad', 'nadp',
             '10thf', 'amet', 'chor', 'coa', 'atp', 'ctp', 'gtp', 'utp', 'enter',
             'gthrd', 'hemeO', 'malcoa', 'mlthf', 'mobd', 'pheme', 'q8', 'mql8',
         'ptrc', 'ribflv', 'sheme', 'spmd', 'succoa', 'udcpdp']

aas.extend(mets)

# Rerun these so that they do not secrete _c -> _e
aas.append('default')

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
anaerobic = False

target_to_flux = {}
target_to_shadow = {}
target_to_reduced = {}
if anaerobic:
    suffix = '_anaerobic'
else:
    suffix = ''
for met in list(aux_to_ko.keys()):

    print(met)
    model = deepcopy(me)
    if anaerobic:
        model.reactions.EX_o2_e.lower_bound = 0

    if met == 'thr__L':
        model.reactions.EX_gly_e.upper_bound = 0


    for r in aux_to_ko[met]:
        for rxn in model.process_data.get_by_id(r).parent_reactions:
            print('knocked out', rxn, 'for', met)
            rxn.knock_out()

    if met == 'default':
        print('Running', met, 'with no modifications')
        me_nlp = ME_NLP1(model, growth_key='mu')
        me_nlp.bisectmu(precision=1e-8, mumax=1.5, mumin=.4996)
        target_to_flux[met] = model.solution.x_dict
        target_to_shadow[met] = dict(
            zip(model.metabolites.list_attr('id'), me_nlp.pi))
        target_to_reduced[met] = dict(
            zip(model.reactions.list_attr('id'), me_nlp.rc))
    else:
        # add exchange for met
        model.add_reaction(cobrame.MEReaction('EX_%s_c' % met))
        model.reactions.get_by_id('EX_%s_c' % met).add_metabolites(
            {'%s_c' % met: -1})
        model.objective = model.reactions.get_by_id('EX_%s_c' % met)

        if met == 'thf':
            r = cobra.Reaction('SK_dhf_c')
            model.add_reaction(r)
            r.add_metabolites({'dhf_c': -1})

        try:
            model.reactions.get_by_id('EX_%s_e' % met).upper_bound = 0
        except:
            print('No extracellular exchange for ', met)

        aux_uptake_r = model.reactions.get_by_id('EX_%s_c' % met)
        aux_uptake_r.lower_bound = -1000

        me_nlp = ME_NLP1(model, growth_key='mu')
        me_nlp.compiled_expressions = me_nlp.compile_expressions()
        x, status, hs = me_nlp.solvelp(.1)

        if status != 'optimal':
            embed()
            print(met, 'does not solve', status)
            target_to_flux[met] = {}
            target_to_shadow[met] = {}
            target_to_reduced[met] = {}
            continue

        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-8, mumax=1.5,
                                                 mumin=0, basis=hs)

        max_uptake = model.solution.x_dict[aux_uptake_r.id]
        titrate_list = np.linspace(max_uptake/20, max_uptake, 20)
        for uptake in titrate_list:
            aux_uptake_r.lower_bound = uptake
            muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-8, mumax=1.5,
                                                     mumin=0, basis=hs)
            col = met + '_' + str(abs(uptake))
            target_to_flux[col] = model.solution.x_dict
            target_to_shadow[col] = dict(
                zip(model.metabolites.list_attr('id'), me_nlp.pi))
            target_to_reduced[col] = dict(
                zip(model.reactions.list_attr('id'), me_nlp.rc))

    save_dir = 'titration%s' % suffix
    pd.DataFrame(target_to_flux).to_csv('./%s/fluxes.csv' % save_dir)
    pd.DataFrame(target_to_shadow).to_csv('./%s/shadow_wrt_target.csv' %
                                          save_dir)
    pd.DataFrame(target_to_reduced).to_csv('./%s/reduced_costs_wrt_target.csv'
                                           % save_dir)
