def delete_non_active_network(original_cmod, which_zeros='flux_zero', zero_tol=1e-9, opt_tol=1e-8):
    """
    Deletes reactions of SBML-model if zero in FBA-solution.
    Checks if the objective function is really not changed due to the deletions.
    :param original_cmod: cbmpy.CBModel
            An FBA should already have been performed
    :param which_zeros: string
            Determines which reactions are deleted. Options
                'flux_zero' deletes all reactions that have zero flux in the optimum
                'FVA_zero' deletes only reactions that have zero flux in all FVA solutions
    :param zero_tol: float
            Determines when a reaction is zero enough
    :param opt_tol: float
            We consider the objective function as not changing if the change is less than this value
    :return cmod: CBMPy-model
            This model was reduced in size by deleting reactions and metabolites that were not active
    """
    delete_reaction = []  # Initializes list for reactions that need to be deleted
    opt_obj = original_cmod.getObjFuncValue()  # Needed for checking if we did not remove too much
    zero_tol = zero_tol * opt_obj
    deletion_success = False
    counter = 0
    N_ITER = 10
    cmod = original_cmod.clone()

    if which_zeros is 'FVA_zero':
        cbm.doFVA(cmod)

    while not deletion_success:  # Runs until we have thrown away only reactions that do not affect the objective
        for rid in cmod.getReactionIds():
            reaction = cmod.getReaction(rid)
            if abs(reaction.value) <= zero_tol:
                if which_zeros is 'flux_zero':
                    delete_reaction.append(rid)
                elif which_zeros is 'FVA_zero':
                    if max(np.abs([reaction.fva_max, reaction.fva_min])) <= zero_tol:
                        delete_reaction.append(rid)

        # delete all reactions in delete_reaction from model
        for rid in delete_reaction:
            print(rid)
            cmod.deleteReactionAndBounds(rid)

        cbm.doFBA(cmod)
        changed_obj = abs(cmod.getObjFuncValue() - opt_obj)
        if changed_obj <= opt_tol:  # Check if objective value didn't change
            deletion_success = True
        else:
            cmod = original_cmod
            zero_tol = zero_tol / 10  # If we have removed to much: Adjust zero_tol and try again

        if counter <= N_ITER:
            counter += 1
        else:
            print("Reaction deletion did not succeed within %d iterations." % N_ITER)

    # Then delete metabolites that are no longer used by any reaction
    stoich_matrix = cmod.N.array
    n_reacs_per_metab = np.count_nonzero(stoich_matrix, axis=1)
    inactive_metabs = [cmod.species[metab_index].id for metab_index in range(stoich_matrix.shape[0]) if
                       n_reacs_per_metab[metab_index] == 0]

    for metab_id in inactive_metabs:
        print('Deleting %s because it is not used in active network.' % metab_id)
        cmod.deleteSpecies(metab_id)

    return cmod


def create_KO_models(original_cmod, n_KOs=3, zero_tol=1e-6, opt_tol=1e-8):
    """
    Deletes reactions of SBML-model if zero in FBA-solution.
    Checks if the objective function is really not changed due to the deletions.
    :param cmod: cbmpy.CBModel
            An FBA should already have been performed
    :param which_zeros: string
            Determines which reactions are deleted. Options
                'flux_zero' deletes all reactions that have zero flux in the optimum
                'FVA_zero' deletes only reactions that have zero flux in all FVA solutions
    :param zero_tol: float
            Determines when a reaction is zero enough
    :param opt_tol: float
            We consider the objective function as not changing if the change is less than this value
    """
    active_reaction_lists = []
    KOs_done = []
    opt_obj = original_cmod.getObjFuncValue()
    obj_id = original_cmod.obj_func.fluxObjectives[0].reaction
    zero_tol = zero_tol * opt_obj

    # Get reactions that are active in first model
    active_reacs = [rid for rid in original_cmod.getReactionIds() if
                    abs(original_cmod.getReaction(rid).value) > zero_tol]
    active_reaction_lists.append(active_reacs)
    KO_candidates = [rid for rid in active_reacs if not rid[:4] == 'R_EX' and not rid == obj_id]

    for KO_ind in range(n_KOs):
        if not len(KO_candidates):
            break

        test_cmod = original_cmod.clone()
        KO_success = False
        counter = 0
        while not KO_success:
            KOs_to_do = KOs_done + [KO_candidates[counter]]
            for KO_to_do in KOs_to_do:
                test_cmod.deleteReactionAndBounds(KO_to_do)

            cbm.doFBA(test_cmod)
            if test_cmod.getObjFuncValue() > 0.01 * opt_obj:
                KO_success = True
            else:
                test_cmod = original_cmod.clone()
                counter += 1

                if counter == len(KO_candidates):
                    break

        if KO_success:
            active_reacs_KO = [rid for rid in test_cmod.getReactionIds() if
                               abs(test_cmod.getReaction(rid).value) > zero_tol]

            active_reaction_lists.append(list(set(active_reaction_lists[-1]) | set(active_reacs_KO)))
            KOs_done = KOs_to_do
            KO_candidates = list(set(KO_candidates).intersection(set(active_reacs_KO)))

    # Now make models with only the active reactions in the list above
    cmod_list = []
    for model_ind in range(1, len(active_reaction_lists)):
        active_rids = active_reaction_lists[model_ind]
        new_model = original_cmod.clone()
        delete_reactions = [rid for rid in new_model.getReactionIds() if rid not in active_rids]
        for rid in delete_reactions:
            new_model.deleteReactionAndBounds(rid)

        cbm.doFBA(new_model)
        new_model.buildStoichMatrix()
        stoich_matrix = new_model.N.array
        n_reacs_per_metab = np.count_nonzero(stoich_matrix, axis=1)
        inactive_metabs = [new_model.species[metab_index].id for metab_index in range(stoich_matrix.shape[0]) if
                           n_reacs_per_metab[metab_index] == 0]

        for metab_id in inactive_metabs:
            print('Deleting %s because it is not used in active network.' % metab_id)
            new_model.deleteSpecies(metab_id)

        new_model.buildStoichMatrix()
        cmod_list.append(new_model)

    return cmod_list