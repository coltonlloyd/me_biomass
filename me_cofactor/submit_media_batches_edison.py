import os
from multime.auxotroph_analysis import load_model
aux_to_ko = {'default': [],
             'pydxn': ['PDX5PS1', 'PDX5PS2'], # PDX5PS in iJO, but unlumped for ME
             'thm': ['THZPSN31'],
             'nac': ['ASPO3', 'ASPO4', 'ASPO5', 'ASPO6'],
             'thf': ['DHFR'],  # actually GCALLD, but seems unlikely
             'met__L': ['HSST'],  # from flexneri 2a
             'pnto__R': ['PANTS'],
             'ribflv': ['RBFSb'],
             'trp__L': ['IGPS'],
             'leu__L': ['IPMD'],  # from DH10b
             'btn': ['ALLTN', 'DBTS'],
             'phe__L': ['PPNDH'],
             'his__L': ['HISTD'],
             'tyr__L': ['PPND'],

             }

me = load_model.load_me_model(json=True)

aerobicity='aerobic'
gene = ''

if aerobicity == 'anaerobic':
    prefix = '_anaerobic'
else:
    prefix = ''

source_dir = os.getcwd() + '/knockout_sims/'
if not os.path.isdir(source_dir):
    os.mkdir(source_dir)

for auxotrophy in aux_to_ko:
    for source in ["C", "P", "S", "N"]:
        for r in me.reactions.query('EX_'):

            if not r.id.startswith('EX_'):
                continue
            if r.lower_bound < 0:
                print('Skipping and not setting bound of %s to 0' % r.id)
                continue
            met = me.metabolites.get_by_id(r.id.replace('EX_', ''))
            if met.elements.get(source, 0) == 0:
                print('No %s in %s' % (source, met.id))
                continue
            media = r.id

            output_file = auxotrophy + '_' + source + '_' + media
            if os.path.exists(os.path.join(source_dir, gene + '%s_sol.json' % prefix)):
                print(gene, 'already solved')
                continue
            os.system("sbatch edison_submit_job %s %s %s %s %s" %
                      (gene, aerobicity, auxotrophy, media, source))
