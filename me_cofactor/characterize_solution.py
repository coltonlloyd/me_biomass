from collections import defaultdict
from sympy import Basic, Symbol
import pandas as pd
import cobra.test
from cobrame.util import dogma

skip_list = [
    'ComplexFormation']  # ['TranslationReaction']# ['SummaryVariable', 'ComplexFormation',


#   'TranscriptionReaction', 'TranslationReaction']
ijo = cobra.test.create_test_model('ecoli')


def compare_cofactor_to_ijo_biomass(model,
                                    currency_met_to_synthesis_rxn=dict(),
                                    solution=None, growth_norm=True):
    # capture all forms of cofactors in iJO1366 BOF
    met_to_met_forms = {'thf_c': ['thf_c', '10fthf_c', 'mlthf_c', '5mthf_c'],
                        'coa_c': ['coa_c', 'accoa_c', 'succoa_c', 'malcoa_c'],
                        'gthox_c': ['gthrd_c', 'gthox_c'],
                        'q8_c': ['q8_c', 'q8h2_c'],
                        '2dmmq8_c': ['2dmmq8_c', '2dmmql8_c'],
                        'mqn8_c': ['mqn8_c', 'mql8_c'],
                        'ribflv_c': ['ribflv_c', 'fad_c', 'fmn_c']}

    if solution:
        model.solution = solution

    biomass_rxn = ijo.reactions.Ec_biomass_iJO1366_WT_53p95M

    growth_rate = 1 if not growth_norm else model.solution.x_dict['biomass_dilution']
    me_demand = {}
    me_demand['ME_model'] = defaultdict(float)
    me_demand['Measured'] = defaultdict(float)
    x_dict = model.solution.x_dict
    for met_id, r_id in currency_met_to_synthesis_rxn.items():
        # some currency mets are synthesized from others. Just account for
        # synthesis of root cofactor
        if not r_id:
            continue
        met = model.metabolites.get_by_id(met_id + '_c')
        for r_obj in model.process_data.get_by_id(r_id).parent_reactions:
            if 'FWD' not in r_obj.id:
                continue

            me_demand['ME_model'][met_id] += (x_dict[r_obj.id] / growth_rate)

        # must make nad to make nadp
        if met_id == 'nad':
            nadp_r = currency_met_to_synthesis_rxn['nadp']
            for r_obj in model.process_data.get_by_id(nadp_r).parent_reactions:
                if 'FWD' not in r_obj.id:
                    continue
                me_demand['ME_model'][met_id] -= (x_dict[r_obj.id] / growth_rate)
    for met in me_demand['ME_model']:
        for other_met in met_to_met_forms.get(met + '_c', [met + '_c']):
            ijo_met = ijo.metabolites.get_by_id(other_met)

            if ijo_met in biomass_rxn.metabolites and biomass_rxn._metabolites[
                ijo_met] < 0:
                me_demand['Measured'][met] += abs(
                    biomass_rxn._metabolites[ijo_met])

    return pd.DataFrame(me_demand)


def compare_to_ijo_biomass(model, kind='amino_acid', solution=None,
                           growth_norm=True):

    if solution:
        model.solution = solution
    biomass_rxn = ijo.reactions.Ec_biomass_iJO1366_WT_53p95M
    me_demand = defaultdict(float)
    x_dict = model.solution.x_dict

    # These are reactions that incorporate metabolites into biomass
    skip_list = ['SummaryVariable', 'ComplexFormation',
                 'TranscriptionReaction', 'TranslationReaction']

    growth_rate = 1 if not growth_norm else model.solution.x_dict['biomass_dilution']
    mu = Symbol('mu')

    if solution:
        model.solution = solution

    if kind == 'amino_acid':
        for r in model.reactions.query('translation_'):
            aa_count = r.translation_data.amino_acid_count
            for aa_letter, aa in dogma.amino_acids.items():
                me_demand[aa] += (aa_count[aa] * r.x / growth_rate)
    elif kind == 'cofactors':
        for d in model.complex_data:
            for mod, num in d.subreactions.items():
                me_demand[mod.replace('mod_', '')] += \
                    (x_dict[d.formation.id] * num / growth_rate)

    else:
        for met_id in biomass_rxn.metabolites:
            met = model.metabolites.get_by_id(met_id.id)
            for r in met.reactions:
                if r.__class__.__name__ not in skip_list:
                    stoich = r._metabolites[met]
                    if isinstance(stoich, Basic):
                        stoich = stoich.subs(mu, growth_rate)
                    me_demand[met_id.id] += (
                            x_dict[r.id] * stoich / growth_rate)

    compare = dict()
    compare['ME_gr_%.2f' % growth_rate] = me_demand.copy()
    compare['Measured'] = {}
    for met in me_demand:

        if met in ijo.metabolites:
            ijo_met = ijo.metabolites.get_by_id(met)
        else:
            continue

        if ijo_met in biomass_rxn.metabolites and \
                biomass_rxn._metabolites[ijo_met] < 0:
            compare['Measured'][met] = \
                abs(biomass_rxn._metabolites[ijo_met])

    return pd.DataFrame.from_dict(compare).dropna(how='any')