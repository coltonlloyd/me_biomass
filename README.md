# Computing and analyzing ME-model predicted biomass compositions

Scripts to reproduce results in "Computation of condition-dependent proteome allocation reveals variability in the macro and micro nutrient requirements for growth"

## Running simulations

`python submit_media_edison.py`
 - Runs ME-model simulations for all possible carbon, nitrogen, phosphorus, and sulfur sources, aerobically and anaerobically
 
 
 `python run_auxotroph_titration`
 - Runs ME-model simulations for amino acid and cofactor auxotroph from excess to 1/20th the optimal amount of the essential metabolite
 
 ## Perform analysis
 
`part_1.ipynb`
  - Recreate Figures 2 through 5
  
 `part_2.ipynb`
 - Recreate Figure 6
