import os
from multime.auxotroph_analysis import load_model

me = load_model.load_me_model(json=True)

aerobicity='anaerobic'

if aerobicity == 'anaerobic':
    prefix = '_anaerobic'
else:
    prefix = ''
for gene_obj in list(me.translation_data):
    gene = gene_obj.id

    source_dir = os.getcwd() + '/knockout_sims/'
    if not os.path.isdir(source_dir):
        os.mkdir(source_dir)
    if os.path.exists(os.path.join(source_dir, gene + '%s_sol.json' % prefix)):
        print(gene, 'already solved')
        continue
    os.system("sbatch edison_submit_job %s %s" % (gene, aerobicity))
