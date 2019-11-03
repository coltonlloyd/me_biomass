from __future__ import print_function, division, absolute_import

from multiprocessing import Pool
import os

from multime import auxotroph_analysis

from multime.auxotroph_analysis import load_model


sim_script_dir = list(auxotroph_analysis.__path__)[0]


def run_pool(function, values, processes=20):
    pool = Pool(processes=processes)
    pool.map(function, values)


def submit_job(gene):
    os.system("python3.6 %s/simulate_knockout.py %s" %
              (sim_script_dir, gene))

me = load_model.load_me_model(json=True)

genes_to_run = []
for gene_obj in list(me.translation_data):
    gene = gene_obj.id

    source_dir = os.getcwd() + '/knockout_sims/'
    if not os.path.isdir(source_dir):
        os.mkdir(source_dir)
    if os.path.exists(os.path.join(source_dir, gene + '_sol.json')):
        print(gene, 'already solved')
        continue

    genes_to_run.append(gene)

print(genes_to_run)
run_pool(submit_job, genes_to_run)
