
def solve_me_model(me, max_mu, precision=1e-6, min_mu=0, using_soplex=True,
                   compiled_expressions=None):
    if using_soplex:
        from cobrame.solve.algorithms import binary_search
        binary_search(me, min_mu=min_mu, max_mu=max_mu, debug=True,
                      mu_accuracy=precision,
                      compiled_expressions=compiled_expressions,
                      solver='gurobi')
    else:
        from qminospy.me1 import ME_NLP1
        # The object containing solveME methods--composite that uses a ME
        # model object
        me_nlp = ME_NLP1(me, growth_key='mu')
        # Use bisection for now (until the NLP formulation is worked out)
        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision,
                                                 mumax=max_mu)
        me.solution.f = me.solution.x_dict['biomass_dilution']


def show_escher_map(me, solution=None):
    import escher
    view = escher.Builder("iJO1366.Central metabolism")
    view.reaction_data = me.get_metabolic_flux(solution=solution)
    return view