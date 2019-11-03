from qminospy.me1 import ME_NLP1
from copy import deepcopy
import pandas as pd
import pickle
import cobra
import cobrame


# Load and update ME-model
with open('/home/sbrg-cjlloyd/multime/multime/iJL1678b_ML_keffs.pickle', 'rb') as f:
    me = pickle.load(f)

me.reactions.get_by_id(
        'PDH_FWD_PYRUVATEDEH-CPLX_mod_mg2_mod_fad_mod_thmpp_mod_lipo').keff = 1500.
me.reactions.get_by_id(
        'PDH_FWD_PYRUVATEDEH-CPLX_mod_mg2_mod_fad_mod_thmpp_mod_lipo').update()

# Assume the F6PA isozyme has the same stoichiometry as the primary isozyme
# TODO CPLX0-201 has a stoichiometry of 10, not 12
me.process_data.get_by_id('EG11905-MONOMER').stoichiometry['protein_b3946'] = 12
me.reactions.get_by_id('formation_EG11905-MONOMER').update()
print(me.reactions.get_by_id('formation_EG11905-MONOMER').reaction)


aas = [i.replace('_c', '') for i in me.process_data.b2020.amino_acid_count.keys()]
mets = ['quln', 'thf', 'thm', 'fad', 'btn', 'hemed', 'fmn', 'nac', 'thmpp', 'nad', 'nadp',
             '10thf', 'amet', 'chor', 'coa', 'atp', 'ctp', 'gtp', 'utp', 'enter',
             'gthrd', 'hemeO', 'malcoa', 'mlthf', 'mobd', 'pheme', 'q8', 'mql8',
         'ptrc', 'ribflv', 'sheme', 'spmd', 'succoa', 'udcpdp']

aas.extend(mets)

# Rerun these so that they do not secrete _c -> _e
aas.append('default')

aux_to_ko = {'thm': ['THZPSN31'],
             'nac': ['ASPO3', 'ASPO4', 'ASPO5', 'ASPO6'],
             'thf': ['DHFR'],  # actually GCALLD, but seems unlikely
             'met__L': ['HSST'],  # from flexneri 2a
             'trp__L': ['IGPS'],
             'leu__L': ['IPMD']  # from DH10b
             }

sim_kind = 'supplemented'
coupled_currency_mets = ['mlthf', '10fthf', 'nadh', 'nadph', 'rbflvrd', 'coa', 'gthrd']
target_to_flux = {}
target_to_shadow = {}
target_to_reduced = {}
for met in aas:

    print(met)
    model = deepcopy(me)
    # Turn off other growth dependent demands of currency mets
    model.reactions.biomass_constituent_demand.lower_bound = 0
    model.reactions.biomass_constituent_demand.upper_bound = 0
    model.reactions.EX_fe2_e.lower_bound = 0
    model.reactions.MOX_REV_CPLX_dummy.upper_bound = 0
    model.reactions.EX_LI_c.lower_bound = 0
    model.reactions.EX_cs_e.lower_bound = 0
    model.reactions.EX_tl_c.lower_bound = 0
    model.reactions.EX_cd2_e.lower_bound = 0
    model.reactions.EX_cbl1_e.lower_bound = -1000
    model.reactions.tl_generation_FWD_CPLX_dummy.knock_out()

    model.process_data.get_by_id(
        '2OXOGLUTARATEDEH-CPLX_mod_mg2_mod_lipo').subreactions[
        'mod_thmpp_c'] = 1.
    model.process_data.get_by_id(
        '2OXOGLUTARATEDEH-CPLX_mod_mg2_mod_lipo').subreactions[
        'mod_fad_c'] = 1.
    model.reactions.get_by_id(
        'formation_2OXOGLUTARATEDEH-CPLX_mod_mg2_mod_lipo').update()

    model.process_data.get_by_id('ACETOLACTSYNIII-CPLX').subreactions[
        'mod_thmpp_c'] = 1.
    model.process_data.get_by_id('ACETOLACTSYNIII-CPLX').subreactions[
        'mod_fad_c'] = 1.
    model.reactions.get_by_id('formation_ACETOLACTSYNIII-CPLX').update()

    for cur_met in coupled_currency_mets:
        for comp in ['_c', '_p']:
            if cur_met + comp not in model.metabolites:
                print(cur_met + comp, 'not in model')
                continue
            met_obj = model.metabolites.get_by_id(cur_met + comp)
            for r in met_obj.reactions:
                stoich = r._metabolites[met_obj]
                if isinstance(r, cobrame.MetabolicReaction) and stoich < 0:
                    r.add_metabolites({met_obj: -abs(stoich) * cobrame.mu/650.}, combine=True)
                    print(r.reaction)

    if met != 'default':
        model.add_reaction(cobrame.MEReaction('EX_%s_c' % met))
        model.reactions.get_by_id('EX_%s_c' % met).add_metabolites(
            {'%s_c' % met: -1})
        model.objective = model.reactions.get_by_id('EX_%s_c' % met)

    if met == 'default':
        print('Running', met, 'with no modifications')
        me_nlp = ME_NLP1(model, growth_key='mu')
        me_nlp.bisectmu(precision=1e-8, mumax=1.5, mumin=.4996)
    elif sim_kind == 'production':
        me_nlp = ME_NLP1(model, growth_key='mu')
        me_nlp.solvelp(.5)
        print(model.solution.x_dict['EX_%s_c' % met])
        model.reactions.get_by_id('EX_%s_c' % met).lower_bound = \
            model.solution.x_dict['EX_%s_c' % met]
        me_nlp.bisectmu(precision=1e-8, mumax=.51, mumin=.49)
    elif sim_kind == 'production_at_5':
        me_nlp = ME_NLP1(model, growth_key='mu')
        me_nlp.solvelp(.5)
        print(model.solution.x_dict['EX_%s_c' % met])
    elif sim_kind == 'supplemented':
        try:
            model.reactions.get_by_id('EX_%s_e' % met).upper_bound = 0
        except:
            print('No extracellular exchange for ', met)
        me_nlp = ME_NLP1(model, growth_key='mu')
        model.reactions.get_by_id('EX_%s_c' % met).lower_bound = -1000
        me_nlp.bisectmu(precision=1e-8, mumax=1.5, mumin=0)

    else:
        raise Exception('Not valid sim_kind/met combination')

    target_to_flux[met] = model.solution.x_dict
    target_to_shadow[met] = dict(
        zip(model.metabolites.list_attr('id'), me_nlp.pi))
    target_to_reduced[met] = dict(
        zip(model.reactions.list_attr('id'), me_nlp.rc))

    pd.DataFrame(target_to_flux).to_csv('./%s/fluxes.csv' % sim_kind)
    pd.DataFrame(target_to_shadow).to_csv('./%s/shadow_wrt_target.csv' % sim_kind)
    pd.DataFrame(target_to_reduced).to_csv(
        './%s/reduced_costs_wrt_target.csv' % sim_kind)

