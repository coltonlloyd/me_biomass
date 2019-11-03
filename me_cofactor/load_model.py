
import pickle
from cobrame.io.json import load_json_me_model
import cobrame
from sympy import Basic
from os.path import dirname, abspath

currency_met_to_synthesis_rxn = {'coa': 'DPCOAK',
                                 'thf': 'DHFS',
                                 # use this reaction because DHFR is coupled to dUMP synthesis
                                 'nadp': 'NADK',
                                 'nad': 'NADS1',
                                 # need two reactions, one creates intermediate for NADK, other couples
                                 'ribflv': 'RBFSb',
                                 'gthox': 'GTHS',
                                 'q8': 'DMQMT',
                                 '2dmmq8': 'DHNAOT4',
                                 # must make 2dmmql8 to produce mql8
                                 'mqn8': 'AMMQLT8',
                                 'fmn': '',
                                 'fad': ''
                                 }

met_to_name = {'2dmmq8': '2-Demethyl-\nmenaquinone 8',
 '2fe2s': '2Fe-2S',
 '4fe4s': '4Fe-4S',
 'adocbl': 'Adenosyl-\ncobalamin',
 'cbl1': 'Cobalamin',
 'ala__L': 'L-Alanine',
 'arg__L': 'L-Arginine',
 'asn__L': 'L-Asparagine',
 'asp__L': 'L-Aspartate',
 'bmocogdp': 'bis-MGD',
 'btn': 'Biotin',
 'coa': 'Coenzyme A',
 'cys__L': 'L-Cysteine',
 'gln__L': 'L-Glutamine',
 'glu__L': 'L-Glutamate',
 'gly': 'Glycine',
 'gthox': 'Glutathione',
 'hemeO': 'Heme O',
 'his__L': 'L-Histidine',
 'ile__L': 'L-Isoleucine',
 'leu__L': 'L-Leucine',
 'met__L': 'L-Methionine',
 'mqn8': 'Menaquinone 8',
 'nad': 'NAD',
 'nadp': 'NADP',
 'phe__L': 'L-Phenylalanine',
 'pheme': 'Protoheme',
 'pro__L': 'L-Proline',
 'pydx5p': "Pyridoxal \n 5'-phosphate",
 'q8': 'Ubiquinone-8',
 'ribflv': 'Riboflavin',
 'ser__L': 'L-Serine',
 'sheme': 'Siroheme',
 'thf': 'Tetrahydrofolate',
 'thmpp': 'Thiamine \n diphosphate',
 'thr__L': 'L-Threonine',
 'trp__L': 'L-Tryptophan',
 'tyr__L': 'L-Tyrosine',
 'val__L': 'L-Valine'}


here = dirname(abspath(__file__))

def load_me_model(json=False):
    # Load and update ME-model
    if not json:
        with open('/home/sbrg-cjlloyd/multime/multime/iJL1678b_ML_keffs.pickle',
                  'rb') as f:
            model = pickle.load(f)
    else:
        model = load_json_me_model('%s/iJL1678b_ML_keffs.json' % here)

    model.reactions.get_by_id(
        'PDH_FWD_PYRUVATEDEH-CPLX_mod_mg2_mod_fad_mod_thmpp_mod_lipo').keff = 1500.
    model.reactions.get_by_id(
        'PDH_FWD_PYRUVATEDEH-CPLX_mod_mg2_mod_fad_mod_thmpp_mod_lipo').update()

    # Assume the F6PA isozyme has the same stoichiometry as the primary isozyme
    # TODO CPLX0-201 has a stoichiometry of 10, not 12
    model.process_data.get_by_id('EG11905-MONOMER').stoichiometry[
        'protein_b3946'] = 12
    model.reactions.get_by_id('formation_EG11905-MONOMER').update()
    print(model.reactions.get_by_id('formation_EG11905-MONOMER').reaction)

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

    for cur_met in currency_met_to_synthesis_rxn.keys():
        mult = 1e4  # USE THIS VALUE AS MULTIPLIER FOR ALL COFACTORS
        print(mult)
        for comp in ['_c', '_p']:
            if cur_met + comp not in model.metabolites:
                print(cur_met + comp, 'not in model')
                continue
            met_obj = model.metabolites.get_by_id(cur_met + comp)
            for r in met_obj.reactions:
                # Double check that reaction is irreversible
                if r.lower_bound < 0 and r.upper_bound > 0:
                    print('Not changing reversible reaction %s, this might be a sign of a problem' % r.id)
                    continue
                stoich = r._metabolites[met_obj]

                # if stoichiometry of cofactor already includes symoblic
                # term (e.g., tRNA charging reactions) then skip
                if isinstance(stoich, Basic):
                    continue

                # only update coefficient if metabolic or translation reaction,
                # this skips complex formation and tRNA charging reactions
                if not isinstance(r, cobrame.MetabolicReaction) and not isinstance(r, cobrame.TranslationReaction):
                    continue

                # update stoichiometry in forward/reverse direction
                if stoich < 0 and r.upper_bound > 0 and r.lower_bound == 0:
                    r.add_metabolites(
                        {met_obj: -abs(stoich) * cobrame.mu / mult},
                        combine=True)
                elif isinstance(r,
                                cobrame.MetabolicReaction) and stoich > 0 and r.lower_bound < 0 and r.upper_bound == 0:
                    r.add_metabolites(
                        {met_obj: abs(stoich) * cobrame.mu / mult},
                        combine=True)
    return model
