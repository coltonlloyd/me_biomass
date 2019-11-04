import os
from me_cofactor import load_model

me = load_model.load_me_model(json=True)

gene = ''
auxotrophy = 'default'
source_dir = os.getcwd() + '/media_sims/'
if not os.path.isdir(source_dir):
    os.mkdir(source_dir)

for aerobicity in ['aerobic', 'anaerobic']:

    for source in ["C", "P", "S", "N"]:
        for r in me.reactions.query('EX_'):

            if not r.id.startswith('EX_'):
                continue
            met = me.metabolites.get_by_id(r.id.replace('EX_', ''))
            if met.elements.get(source, 0) == 0:
                print('No %s in %s' % (source, met.id))
                continue
            media = r.id

            output_file = aerobicity + '_' + source + '_' + media
            if os.path.exists(os.path.join(source_dir, '%s_sol.json' % output_file)):
                print(output_file, 'already solved')
                continue
            os.system("sbatch edison_submit_job %s %s %s %s %s" %
                      (gene, aerobicity, auxotrophy, media, source))