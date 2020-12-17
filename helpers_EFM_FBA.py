import copy
import os
from fractions import Fraction

import cbmpy as cbm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from matplotlib import cm
from scipy.optimize import linprog

from ecmtool import extract_sbml_stoichiometry, get_conversion_cone
from ecmtool.helpers import to_fractions, unique
from ecmtool.network import add_reaction_tags  # was updated in v.0.1.6, was add_debug_tags


def findConstraintMetabolites(nzrc_dictionaries, intermediate_cmod):
    """
    Checks which metabolite is coupled to the constrained reactions. If no metab is coupled, it stores the reaction
    in reactions_to_tag, indicating that a virtual metabolite should be added to track this constrained reaction
    using conversions
    :param nzrc_dictionaries: list of dictionaries
            Dictionary with 4 entries: (reaction id, reduced cost, obj/sec_obj/cons, value of bound that was hit)
    :param intermediate_cmod: CBMPy-model object
            metabolic model
    :return nzrc_dictionaries: list of dictionaries
            We have added to the dictionaries 5) species-ID of metabolite that is coupled to the constrained reaction
            6) Stoichiometric coefficient with which this species is coupled to the constrained reaction
    :return reactions_to_tag: list of reaction-IDs
    """
    reactions_to_tag = []

    for nzrc_ind, nzrc_dict in enumerate(nzrc_dictionaries):
        re = intermediate_cmod.getReaction(nzrc_dict['rid'])  # Get corresponding reaction
        re_stoich = re.getStoichiometry()

        # Try to find external species used in the reaction
        found_external = False
        for stoich_pair in re_stoich:
            species = intermediate_cmod.getSpecies(stoich_pair[1])
            if species.compartment in ['e', 'Extra_organism']:
                nzrc_dictionaries[nzrc_ind]['ext_metab'] = stoich_pair[1]
                nzrc_dictionaries[nzrc_ind]['ext_metab_stoich'] = stoich_pair[0]
                found_external = True
                break
                # If there are two externals in the reaction, we only need to track one of them

        # If not we add the reaction to the list reactions_to_tag. The name of the virtual metabolite is now based on
        # what is used in ECMtool. I don't know if this can be done more generically
        # TODO: Try to find a more generic way
        if not found_external:
            print('No external is found. This reaction is tagged with a virtual external metabolite')
            rea_id = nzrc_dict['rid']
            reactions_to_tag.append(rea_id)
            nzrc_dictionaries[nzrc_ind]['ext_metab'] = 'virtual_tag_%s' % rea_id
            nzrc_dictionaries[nzrc_ind]['ext_metab_stoich'] = 1

    return nzrc_dictionaries, reactions_to_tag


def get_info_objectives_constraints(nzrc_dictionaries, intermediate_cmod):
    """

    :param nzrc_dictionaries: list of dictionaries
            Dictionary with 4 entries: (reaction id, reduced cost, obj/sec_obj/cons, value of bound that was hit)
            5) species-ID of metabolite that is coupled to the constrained reaction 6) Stoichiometric coefficient
            with which this species is coupled to the constrained reaction
    :param intermediate_cmod: CBMPy model object
    :return infos_obj: list of dictionaries
            Contains the same info as nzrc_dictionaries but only about the objective and the secondary objectives
    :return infos_cons: list of dictionaries
            Contains the same info as nzrc_dictionaries but only about the constraints
    """
    infos_obj = []
    infos_cons = []

    active_obj_rid = intermediate_cmod.getActiveObjectiveReactionIds()[0]
    # TODO: Make the objective metabolite name: 'objective' more generic. Now this is based on that I know what ECMtool adds

    # First create dictionary with information about the objective function.
    obj_dict = {'rid': active_obj_rid, 'nzrc': np.nan, 'obj_cons': 'obj', 'flux_val': np.nan, 'ext_metab': 'objective',
                'ext_metab_stoich': 1.0}
    infos_obj.append(obj_dict)

    # Then separate the dictionaries in nzrc_dictionaries into objective dictionaries and constrained dictionaries
    for i, nzrc_dict in enumerate(nzrc_dictionaries):
        infos_cons.append(nzrc_dict) if nzrc_dict['obj_cons'] == 'cons' else infos_obj.append(nzrc_dict)

    return infos_obj, infos_cons


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
    if which_zeros is 'FVA_zero':
        cbm.doFVA(original_cmod)
    cmod = original_cmod.clone()

    # if which_zeros is 'FVA_zero':
    #     cbm.doFVA(cmod)

    while not deletion_success:  # Runs until we have thrown away only reactions that do not affect the objective
        print("This is round", counter, "of", N_ITER)
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
            cmod.deleteReactionAndBounds(rid)

        cbm.doFBA(cmod)
        changed_obj = abs(cmod.getObjFuncValue() - opt_obj)
        if changed_obj <= opt_tol:  # Check if objective value didn't change
            deletion_success = True
            print("Reaction deletion succeeded.")
        else:
            cmod = original_cmod.clone()
            zero_tol = zero_tol / 10  # If we have removed to much: Adjust zero_tol and try again
            delete_reaction = []

        if counter <= N_ITER:
            counter += 1
        else:
            print("Reaction deletion did not succeed within %d iterations." % N_ITER)
            break

    if deletion_success:
        # Then delete metabolites that are no longer used by any reaction
        stoich_matrix = cmod.N.array
        n_reacs_per_metab = np.count_nonzero(stoich_matrix, axis=1)
        inactive_metabs = [cmod.species[metab_index].id for metab_index in range(stoich_matrix.shape[0]) if
                           n_reacs_per_metab[metab_index] == 0]

        for metab_id in inactive_metabs:
            print('Deleting %s because it is not used in active network.' % metab_id)
            cmod.deleteSpecies(metab_id)

        return cmod
    else:
        return None


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
    print(opt_obj)
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
        if new_model.SOLUTION_STATUS == 'LPS_OPT':
            new_model.buildStoichMatrix()
            stoich_matrix = new_model.N.array
            n_reacs_per_metab = np.count_nonzero(stoich_matrix, axis=1)
            inactive_metabs = [new_model.species[metab_index].id for metab_index in range(stoich_matrix.shape[0]) if
                               n_reacs_per_metab[metab_index] == 0]

            for metab_id in inactive_metabs:
                print('Deleting %s because it is not used in active network.' % metab_id)
                new_model.deleteSpecies(metab_id)

            new_model.buildStoichMatrix()
            print('objective value of model is %d' % new_model.getObjFuncValue())
            cmod_list.append(new_model)
        else:
            print("FBA of this model didn't result in optimal solution; model is skipped.")

    return cmod_list


def get_nzrc(cmod):
    """
    Gets all reactions with non-zero-reduced-costs (nzrc). A reaction has nzrc if the objective function value changes
    if the bounds on this constraint changes. We will only consider reactions with nzrc and flux value not equal to 0.
    These are irreversible reactions that would be favorable as reversible, but we don't want to consider those.
    We also round the nzrc values to 10 decimals, this prevents very small nzrc to pop-up, who won't be present in the
    active network or active FVA network.
    We note if the nzrc is due to a fixed objective flux (e.g. ATP-maintenance) or due to an actual constraint
    :returns nzrc_dictionaries: list of dictionaries
            Dictionary with 4 entries: (reaction id, reduced cost, obj/sec_obj/cons, value of bound that was hit)
    :returns n_objectives: int
            Shows how many objectives there are in the problem. Certainly 1 (objective function), but can be more if
            there are fluxes that should have some value. Reducing these values would increase the objective flux.
    :param cmod: cbmpy.CBModel
            An FBA should already have been performed
    """
    nzrc_dictionaries = []

    # We only consider reactions with nzrcs with flux value not equal to 0.
    # The following finds tuples of reaction ids and their reduced costs.
    non_zero_reduced_cost_pairs = [(rid, cmod.getReaction(rid).reduced_cost) for rid in
                                   cmod.getReactionIds() if
                                   round(abs(cmod.getReaction(rid).reduced_cost), 15) > 0]  # round to 10 decimals

    # TODO: This does not work with more general constraints than flux bound yet.
    n_objectives = 1
    # We assume that we at least have one objective, but lower bounds on positive fluxes can be seen as secondary
    # objectives, since the model needs to satisfy this bound as cheap as possible
    for nzrc_pair in non_zero_reduced_cost_pairs:
        rid = nzrc_pair[0]
        nzrc = nzrc_pair[1]
        flux_value = cmod.getReaction(rid).getValue()
        if not flux_value == 0.:  # We don't consider reactions that are zero.
            if nzrc_pair[1] * flux_value <= 0:  # An extra optimisation: to make this flux as cheap as possible
                n_objectives += 1
                obj_cons = 'sec_obj'
            else:  # A constraint
                obj_cons = 'cons'

            nzrc_dict = {'rid': rid, 'nzrc': nzrc, 'obj_cons': obj_cons, 'flux_val': flux_value}
            nzrc_dictionaries.append(nzrc_dict)

    return nzrc_dictionaries, n_objectives


def get_network_class(file_path, reactions_to_tag=[], use_external_compartment=None):
    """
    :return network: network class
    :param file_path: string
            String with path to the file.
    :param reactions_to_tag: list of strings
            Strings with reaction IDs that should be tagged with a virtual metabolite
    """
    # Stap 1: Function from ECMtool that builds network class
    network = extract_sbml_stoichiometry(file_path, determine_inputs_outputs=True,
                                         use_external_compartment=use_external_compartment)

    # Find indices of reactions that should be tagged
    indices_to_tag = []
    if len(reactions_to_tag) > 0:
        for rid in reactions_to_tag:
            ind_to_tag = [ind for ind, reaction in enumerate(network.reactions) if reaction.id == rid]
            indices_to_tag.append(ind_to_tag[0])

        add_reaction_tags(network, reactions=indices_to_tag)  # Function from ECMtool that adds the virtual metabolites

    return network


def calc_ECMs(file_path, reactions_to_tag=[], print_results=False, hide_metabs=[], use_external_compartment=None,
              only_rays=False):
    """
    Calculates ECMs using ECMtool
    :return ecms: np.array
            This array contains the ECMs as columns and the metabolites as rows
    :param file_path: string
            String with path to the SBML-file.
    :param reactions_to_tag: list with strings
            List with reaction-IDs of reactions that need to be tagged
    :param print_results: Boolean
    :param use_external_compartment=None default # if a string is given, this is used to detect external metabolites
    :param hide_metabs: indices of metabolites that should be ignored
    """
    # Stap 1: netwerk bouwen
    network = extract_sbml_stoichiometry(file_path, determine_inputs_outputs=True,
                                         use_external_compartment=use_external_compartment)
    indices_to_tag = []
    # print([reaction.id for reaction in network.reactions])
    if len(reactions_to_tag) > 0:
        for rid in reactions_to_tag:
            ind_to_tag = [ind for ind, reaction in enumerate(network.reactions) if reaction.id == rid]
            # print(rid)
            # print(ind_to_tag)
            # indices_to_tag.append(ind_to_tag[0])
            indices_to_tag += ind_to_tag

        add_reaction_tags(network, reactions=indices_to_tag)

    if len(hide_metabs) > 0:
        network.hide(hide_metabs)

    full_network = copy.deepcopy(network)
    orig_N = network.N

    # Stap 2: compress network
    if print_results:
        print("\n", "Compressing network")
    network.compress(verbose=True, cycle_removal=False)  # verbose was True

    # Stap 3: Ecms enumereren
    # TODO: add timer on enumerating ECMs. tic toc?
    if print_results:
        print("\n", "Enumerating ECMs")
    cone = network.uncompress(
        get_conversion_cone(network.N, network.external_metabolite_indices(), network.reversible_reaction_indices(),
                            network.input_metabolite_indices(), network.output_metabolite_indices(), only_rays=only_rays,
                            verbose=True))  # verbose was True

    if print_results:
        print_ECMs(cone, indices_to_tag, network, orig_N, add_objective_metabolite=True, check_feasibility=True)

    cone = cone.transpose()  # columns will be the different ECMs, rows are metabolites

    return cone, full_network


def normalize_to_row(matrix, row_ind, not_normalized_yet):
    """
    :param matrix: np.array
            Matrix that should be normalized
    :param row_ind: int
            Row that should be normalized to
    :param not_normalized_yet: list of ints
            List of column indices that still need normalization
    :return: matrix: np.array
            Normalized matrix
    :return: not_normalized_yet: list of ints
            Updated list of column indices that still need normalization
    """
    obj_row = matrix[row_ind, :]
    div_factors = [Fraction(1, 1)] * matrix.shape[1]  # By default, divide by 1
    normalized_indices = []
    ecms_to_be_normalized = []

    # Find those colunns that are not yet normalized, but can be normalized using this row
    for col_ind, ecm_ind in enumerate(not_normalized_yet):
        if obj_row[ecm_ind] != 0:
            div_factors[ecm_ind] = abs(obj_row[ecm_ind])  # If column can be normalized, divide by the obj_row-value
            normalized_indices.append(col_ind)
            ecms_to_be_normalized.append(ecm_ind)

    div_factors = np.array(div_factors)[ecms_to_be_normalized]
    divisor_matrix = np.tile(div_factors, (matrix.shape[0], 1))
    matrix[:, ecms_to_be_normalized] = np.divide(matrix[:, ecms_to_be_normalized], divisor_matrix)
    not_normalized_yet = np.delete(not_normalized_yet, normalized_indices)
    return matrix, not_normalized_yet


# def normalize_to_row_splitted(matrix, row_ind, not_normalized_yet):
#     """
#     Normalizes to row, but also splits a large matrix to circumvent memory issues.
#     :param matrix: np.array
#             Matrix that should be normalized
#     :param row_ind: int
#             Row that should be normalized to
#     :param not_normalized_yet: list of ints
#             List of column indices that still need normalization
#     :return: matrix: np.array
#             Normalized matrix
#     :return: not_normalized_yet: list of ints
#             Updated list of column indices that still need normalization
#     """
#     obj_row = matrix[row_ind, :]
#     div_factors = [Fraction(1, 1)] * matrix.shape[1]  # By default, divide by 1
#     normalized_indices = []
#
#     # Find those colunns that are not yet normalized, but can be normalized using this row
#     for col_ind, ecm_ind in enumerate(not_normalized_yet):
#         if obj_row[ecm_ind] != 0:
#             div_factors[ecm_ind] = abs(obj_row[ecm_ind])  # If column can be normalized, divide by the obj_row-value
#             # eunice edit: was without abs()
#             normalized_indices.append(col_ind)
#
#     not_normalized_yet = np.delete(not_normalized_yet, normalized_indices)
#
#     divisor_matrix = np.tile(div_factors, (matrix.shape[0], 1))
#
#     matrix_old = matrix.copy()
#
#     # TODO: split matrix in sets of 100 (?) columns, now we run into memory issues with large sets of ECMs.
#     for counter in range((matrix.shape[1]- matrix.shape[1]%100)/100):
#         # subset matrix
#         print(list(range(counter*100, (counter+1)*100)))
#         matrix_old
#
#         # np.divide(matrix, divisor_matrix)
#
#
#         # combine processed matrices again by adding them in the same order?
#
#
#     # do the rest of the matrix.shape[1]%100 columns
#     # combine the rest with the big part.
#     matrix = np.divide(matrix, divisor_matrix)
#
#     return matrix, not_normalized_yet

def igen(a, n, m):
    """
    Generates a generator object that splits up a np.array in smaller np.arrays of size n x m, with the remaining in
    a smaller np.array.
    :param a:   np.array
                matrix that will be splitted in smaller matrices
    :param n:   int
                size of matrix length
    :param m    int
                size of matrix width
    """
    i_ = np.arange(a.shape[0]) // n
    j_ = np.arange(a.shape[1]) // m
    for i, j in product(np.unique(i_), np.unique(j_)):
        yield (i, j), a[i_ == i][:, j_ == j]


def matrix_splitter(matrix, split_num):
    col_divider = list(range(matrix.shape[1]))
    dict_col_divider = dict()

    for i in range(int(np.ceil(matrix.shape[1] / split_num))):
        dict_col_divider[i] = col_divider[i * split_num: (i + 1) * split_num]

    dict_matrix = dict()
    for i in dict_col_divider.keys():
        dict_matrix[i] = matrix[:, dict_col_divider[i]]

    return dict_matrix


def normalize_ECMS(ecms_matrix, network, normalization_order=[]):
    """
    Normalizes ECMs first to first objective. If there are also lower bounds that act as a kind of second objective.
    Then normalize the ECMs with zero objective to this second objective.
    :return ecms_matrix: np.array
            This array contains the normalized ECMs as columns and the metabolites as rows
    :param ecms_matrix:
            This array contains the ECMs as columns and the metabolites as rows
    :param network: Network class
            Network class as comes from ECMtool
    :param normalization_order: list of metab-IDs
            List of ordered metabolite-IDs
    """
    MATRIX_SPLIT = 100000
    # Determine an order for normalizing the metabolites, if none is given
    # ECMs that are used often are picked first for normalization. To compare two ECM results, make sure to pick the
    # same normalization order
    if not len(normalization_order):
        normalization_order = determine_normalization_order(ecms_matrix, network)

    not_normalized_yet = list(range(ecms_matrix.shape[1]))  # This tracks which ECMs need to be normalized still
    zero_cols = np.where(np.count_nonzero(ecms_matrix, axis=0) == 0)[0]
    not_normalized_yet = np.delete(not_normalized_yet, zero_cols)

    # If matrix.shape[1] > 100000: split the matrix in smaller subsets.
    if ecms_matrix.shape[1] > MATRIX_SPLIT:  # 100000
        # split matrix in smaller subsets
        print('Split ECMs matrix, because number of ECMs is higher then '  + str(MATRIX_SPLIT) + ': ' + str(ecms_matrix.shape[1]))
        print('This takes time.')
        # matrix_dict = dict(igen(ecms_matrix, ecms_matrix.shape[0], 100000))
        matrix_dict = matrix_splitter(ecms_matrix, MATRIX_SPLIT)
        # not_normalized_yet_dict =
        print('ECMs matrix is splitted in ' + str(len(matrix_dict)))

        # Create dict to keep track of which ECMs need to be normalized still.
        dict_not_normalized_yet = dict()
        for i in matrix_dict.keys():
            dict_not_normalized_yet[i] = list(range(matrix_dict[i].shape[1]))

        #  Normalize ECMs
        for key in matrix_dict.keys():
            print('Normalizing matrix ' + str(key) + ' of ' + str(len(matrix_dict)))

            # Then normalize all ECMs to one of the metabolites
            for metab in normalization_order:
                print('Normalizing with ' + metab + ' because ' + str(
                    len(dict_not_normalized_yet[key])) + ' ecms are not yet normalized.')
                if not len(dict_not_normalized_yet[key]):
                    break

                metab_index = [index for index, met in enumerate(network.metabolites) if met.id == metab][0]
                # print(matrix_dict[key])
                matrix_dict[key], dict_not_normalized_yet[key] = normalize_to_row(matrix_dict[key], metab_index,
                                                                                  dict_not_normalized_yet[key])
                # print(matrix_dict[key])
                if len(dict_not_normalized_yet[key]) == 0:
                    break

        # If all ECMs are normalized, i.e. all lists in dict_not_normalized_yet are empty,
        # then normalization is finished; normalized ECMs matrix will be returned
        print('Normalization is done, splitted matrices will be concatenated and returned.')
        # Combine splitted matrix into one again.
        normalized_ecms_matrix = np.concatenate(
            [matrix_dict[i] for i in list(range(len(matrix_dict)))],
            axis=1)
        return normalized_ecms_matrix

    else:
        # Then normalize all ECMs to one of the metabolites
        for metab in normalization_order:
            print('Normalizing with ' + metab + ' because ' + str(
                len(not_normalized_yet)) + ' ecms are not yet normalized.')
            if not len(not_normalized_yet):
                break

            metab_index = [index for index, met in enumerate(network.metabolites) if met.id == metab][0]
            ecms_matrix, not_normalized_yet = normalize_to_row(ecms_matrix, metab_index, not_normalized_yet)

            # If all ECMs are normalized, we can stop normalizing and normalized ecms_matrix will be returned
            if len(not_normalized_yet) == 0:
                return ecms_matrix


def normalize_ECMS_objective_first(ecms_matrix, network, infos_obj, verbose=True):
    """
    Normalizes ECMs first to first objective. If there are also lower bounds that act as a kind of second objective.
    Then normalize the ECMs with zero objective to this second objective.
    :return ecms_matrix: np.array
            This array contains the normalized ECMs as columns and the metabolites as rows
    :param ecms_matrix:
            This array contains the ECMs as columns and the metabolites as rows
    :param network: Network class
            Network class as comes from ECMtool
    :param infos_obj: list of dictionaries
            Dictionaries with information about the different objective functions
    """
    # We first normalize for the objective metabolite
    objective_metab = [info_dict['ext_metab'] for info_dict in infos_obj if np.isnan(info_dict['flux_val'])][0]
    # Then for the metabolites that correspond to another objective (often a lower bound on a certain flux)
    secondary_objectives = [info_dict['ext_metab'] for info_dict in infos_obj if
                            info_dict['ext_metab'] is not objective_metab]
    # Then for all other metabolites
    tertiary_objectives = []
    tertiary_objectives_usage = []
    for metab_ind, metab in enumerate(network.metabolites):
        if metab.id is not objective_metab and metab.id not in secondary_objectives:
            tertiary_objectives.append(metab.id)  # Store all other metabolite ids
            tertiary_objectives_usage.append(
                np.count_nonzero(ecms_matrix[metab_ind, :]))  # And their number of occurrences in ECMs

    # Order the tertiary objectives for how often they are used
    tertiary_objectives = [x for (y, x) in
                           sorted(zip(tertiary_objectives_usage, tertiary_objectives), key=lambda pair: pair[0],
                                  reverse=True)]
    if verbose:
        print('Normalizing in the following order:')
        print(objective_metab)
        print(*secondary_objectives, sep=', ')  # Eunice edit
        print(*tertiary_objectives, sep=', ')  # Eunice edit

    # First normalize everything that can be normalized for the real objective
    ecms_matrix = normalize_ECMS(ecms_matrix, network, [objective_metab] + secondary_objectives + tertiary_objectives)

    return ecms_matrix


def calc_EFMs(network, result_dir, verbose=True):
    """
    Calculates EFMs using EFMtool in Matlab. Saves a csv-file with the EFMs as columns.
    IMPORTANT: This function will only work if the steps on the following site are used:
    https://nl.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
    :return efms_df: Pandas dataframe
            This dataframe contains the EFMs as columns and the reactions as rows
    :param network: Network class
            Network-object as used in ECMtool for which EFMs should be calculated
    """
    import matlab.engine

    # Add your own folder here where efmtool is located
    # efmtool_folder = "C:\\Users\\Eunice van Pelt\\Desktop\\whole-cell-model\\eunice\\Daan\\efmtool"
    efmtool_folder = "C:\\Users\\Daan\\surfdrive\\PhD\\Software\\efmtool"

    # Use the same stoichiometric matrix as used for the ECM-calculation
    stoich_matrix = network.N
    n_orig_reacs = len(network.reactions)
    # Find reaction names in the order used in the network class
    reac_names = [reac.id for reac in network.reactions]
    # Find the same reversibilities as used for the ECM-calculation
    reversibilities = [1 if reac.reversible else 0 for reac in network.reactions]

    # Add an exchange reaction for external metabolites, otherwise no steady-state can be obtained
    for ind_metab, metab in enumerate(network.metabolites):
        if metab.is_external:
            # Determine reversibilities of the exchange reaction, based on information if the metab
            # is an input, output or both
            reversible = 1 if metab.direction == 'both' else 0
            # Create new column (reaction) in the stoichiometric matrix
            col = to_fractions(np.zeros(shape=(stoich_matrix.shape[0], 1)))
            col[ind_metab, 0] = -1 if metab.direction == 'output' else 1
            stoich_matrix = np.append(stoich_matrix, col, axis=1)
            reversibilities.append(reversible)
        #     if metab.is_external:
        #         reversible = 0
        #         if metab.direction in ['input','both']:
        #             col = to_fractions(np.zeros(shape=(stoich_matrix.shape[0], 1)))
        #             col[ind_metab, 0] = 1
        #             stoich_matrix = np.append(stoich_matrix, col, axis=1)
        #             reversibilities.append(reversible)
        #         if metab.direction in ['output', 'both']:
        #             col = to_fractions(np.zeros(shape=(stoich_matrix.shape[0], 1)))
        #             col[ind_metab, 0] = -1
        #             stoich_matrix = np.append(stoich_matrix, col, axis=1)
        #             reversibilities.append(reversible)

    stoich_matrix = stoich_matrix.astype(dtype=float)  # Matlab doesn't deal with the python Fractions

    # Save the stoichiometric matrix and reversibilities in the result directory, so that Matlab can use it
    np.savetxt(os.path.join(result_dir, "stoich_matrix.csv"), stoich_matrix, delimiter=",")
    np.savetxt(os.path.join(result_dir, "reversibilities.csv"), reversibilities, delimiter=",")

    engine = matlab.engine.start_matlab()
    # Uses our own Matlab-function to set EFMtool to the task of calculating EFMs
    result = engine.calculate_efms_wo_SBML(efmtool_folder, result_dir)
    # print("matlab result \n", result)
    size_efms = result.size
    # print(size_efms)
    # The following is a fast way of importing large arrays from Matlab
    efms_matrix = np.reshape(np.array(result._data), size_efms, order='F')
    # Crop virtual exchange reactions of
    efms_matrix = efms_matrix[0:n_orig_reacs, :]
    # Create dataframe with EFMs as the columns, and reactions_names as the rows
    efms_df = pd.DataFrame(efms_matrix, index=reac_names)
    engine.quit()

    return efms_df


def plot_costs(model_dict, infos_obj, infos_cons, cons_IDs=[], obj_val=0.33, show_active=True, result_dir=None, vector_focus=True):
    """
    Plots the costs of all ECMs in cost_df in the cost vector figure, and shows ECM_usage in the FBA-solution by having
    vectors.
    :param result_dir: directory path
            Directory for storing the figures
    :param model_dict: dictionary
            Dictionary with all information about one of the models that we are considering
    :param infos_obj: list
            List of dictionaries concerning the different objectives
    :param infos_cons: list
            List of dictionaries concerning the different constraints
    :param cons_IDs: list of cons_IDs
            List of strings reflecting the different constraints that should be plotted
    :param obj_val: float
            Objective value needed for potentially rescaling the plot
    :param show_active: Boolean
            Show or not show ECM usage
    :param vector_focus: Boolean
            Focus on vectors or show all ECMs
    """
    # Some important numbers
    n_objectives = len(infos_obj)
    scatter_active = 100
    alpha_active = 0.5
    scatter_normal = 15
    quiverwidth = 0.01

    # If the objective is really large, we will show the costs to produce more than one unit objective. We use the
    # scaling factor to calculate the new costs
    scale_obj_val = np.floor(np.log2(obj_val))
    scaling_factor = 2 ** scale_obj_val

    # Find the indices of the ECMs that are active in the FBA solution
    actECMs_inds = list(model_dict['table_cons_df'].loc[model_dict['table_cons_df']['active'] > 0]['orig_ECM_number'])

    # We use gray for normal non-active ECMs, and a different colour for each active ECM
    cmap = cm.get_cmap('tab10', len(actECMs_inds))

    # Make sure that we always have two constraints for the plot
    cons_IDs = select_cons_IDs(infos_cons, cons_IDs)

    # Create one figure with n subplots, where n is the number of 'objectives'. A lower bound that should be met is also
    # considered an objective here.
    if n_objectives > 3:
        fig, axes = plt.subplots(np.int(np.ceil(n_objectives/3)), 3, figsize=(7,6))
        quiverwidth = quiverwidth * 3
    else:
        fig, axes = plt.subplots(1, n_objectives, figsize=(7, 4))
    if n_objectives == 1:
        axes = [axes]  # Necessary evil
    # We loop over the subplots, i.e., the number of objectives
    for index, objective_dict in enumerate(infos_obj):
        if n_objectives > 3:
            ax = axes[np.int(index/3)][index%3]
        else:
            ax = axes[index]
        plt.sca(ax)
        obj_string = objective_dict['ext_metab']

        # In the x_ulims and y_ulims we will store a number of candidates for the upper bound of the plot. We will
        # eventually take the maximum of these candidates
        x_ulims = []
        y_ulims = []
        x_llims = []
        y_llims = []

        # Make cost table with only the constraints that are to be shown
        cost_df, cons_indices = get_cost_df(model_dict, infos_cons, cons_IDs, index, objective_dict,
                                            scaling_factor=scaling_factor)

        if show_active and 'active' not in cost_df.columns:
            print("Can't show active EFMs if this information is not provided.")
            return

        x_label = cons_IDs[0]
        y_label = cons_IDs[1]

        # Determine the bounds for the axes along which we want to plot the costs. It is a bit complicated, but
        # works
        if objective_dict['obj_cons'] == 'obj':
            x_llims.append(-0.05)
            x_ulims.append(max(1.1, min(5 * scaling_factor / obj_val, max(cost_df.values[:, 0]) * 1.1)))
            if not vector_focus:
                x_ulims.append(max(cost_df[cons_IDs[0]]) * 1.1)
                x_llims.append(min(0, min(cost_df[cons_IDs[0]])))
            y_llims.append(-0.05)
            y_ulims.append(max(1.1, min(5 * scaling_factor / obj_val, max(cost_df.values[:, 1]) * 1.1)))
            if not vector_focus:
                y_ulims.append(max(cost_df[cons_IDs[1]]) * 1.1)
                y_llims.append(min(0, min(cost_df[cons_IDs[1]])))
        else:
            x_ulims.append(1.1)
            y_ulims.append(1.1)
            x_llims.append(-0.05)
            if not vector_focus:
                x_ulims.append(max(cost_df[cons_IDs[0]]) * 1.1)
                x_llims.append(min(0, min(cost_df[cons_IDs[0]])))
            y_llims.append(-0.05)
            if not vector_focus:
                y_ulims.append(max(cost_df[cons_IDs[1]]) * 1.1)
                y_llims.append(min(0, min(cost_df[cons_IDs[1]])))

        if 'virtual' in cons_IDs:
            y_ulims = [0.3]
            y_llims = [-0.3]

        # Plot cost dots of normal ECMs
        cost_df.plot.scatter(x=x_label, y=y_label, color='grey', ax=ax,
                             s=scatter_normal, label=model_dict['model_name']) #, legend=False)

        # Plot cost dots of active ECMs
        actECMs = cost_df.loc[cost_df['orig_ECM_number'].isin(actECMs_inds)]
        # The number of active ECMs for this objective can be lower than the number of active ECMs in total, because
        # we cannot plot the costs for an ECM that does not contribute anything to this objective.
        n_actECMs_curr_obj = len(actECMs)

        # Here we will store the x,y-length of the cost vectors if multiplied by its activity
        ecm_usage = np.zeros((n_actECMs_curr_obj, 2))
        # x-length
        ecm_usage[:, 0] = actECMs.values[:, cons_indices[0]] * actECMs['active']
        # y-length
        ecm_usage[:, 1] = actECMs.values[:, cons_indices[1]] * actECMs['active']

        # Construct two columns of starting positions for the vectors for the vector addition with which we show the
        # usage of the ECMs. The end point of the first vector is the starting point of the second vector, etc.
        quiver_start = ecm_usage.cumsum(axis=0)
        # Append starting point (0,0) of the first vector
        quiver_start = np.append(np.zeros((1, 2)), quiver_start, axis=0)
        # Delete last row because this is actually the end point of the last vector
        quiver_start = quiver_start[:-1, :]
        # Add the x,y-lengths of the vectors as two additional columns
        ecm_usage = np.append(quiver_start, ecm_usage, axis=1)

        # The following thus are the x,y-begin positions of the vectors, and the x,y-length
        Xusage, Yusage, Uusage, Vusage = zip(*ecm_usage)

        # Get the indices of the ECMs that contribute to this objective, needed for getting the right colour
        curr_act_inds = [counter for counter, ind in enumerate(actECMs_inds) if ind in actECMs['orig_ECM_number']]
        cmap_curr = cmap.colors[curr_act_inds, :]

        # Plot the costs per unit objective of the active ECMs
        actECMs.plot(kind='scatter', x=x_label, y=y_label, color=cmap_curr, ax=ax,
                     s=scatter_active, alpha=alpha_active, label='active ECM') #, legend=True, label=[str(i) for i in actECMs_inds])

        # Plot the usage of ECMs as vectors
        # for i in range(len(Xusage)):
        #    color = cmap_curr[i, :]
        #    ax.quiver(Xusage[i], Yusage[i], Uusage[i], Vusage[i], pivot='tail', angles='xy', scale_units='xy',
        #              linestyle='--', color=color, scale=1, width=quiverwidth)

        # Set dimensions of the plot and configure axes
        # Also we make a light grey "constraint-box". The allowed region for solutions should be in this box.
        if vector_focus:
            ax.set(adjustable='box', aspect='equal')
        y_ulim = max(y_ulims)
        x_ulim = max(x_ulims)
        x_llim = min(x_llims)
        y_llim = min(y_llims)
        ax.plot([1, 1], [y_llim, y_ulim], '--', color='grey', linewidth=2.0)
        ax.set_xlim([x_llim, x_ulim])
        ax.set_ylim([y_llim, y_ulim])
        # plt.xlabel('Needed fraction of constraint: ' + x_label)
        ax.set(xlabel='')
        if 'virtual' in cons_IDs:
            ax.axes.get_yaxis().set_visible(False)
            ax.axhline(y=0, color='grey', linewidth=2.0)
        else:
            ax.set(ylabel='')
            # plt.ylabel('Needed fraction of constraint: ' + y_label)
            ax.plot([0, 1], [0, 1], '-.', color='grey')
            ax.axes.get_yaxis().set_visible(True)
            ax.plot([x_llim, x_ulim], [1, 1], '--', color='grey', linewidth=2.0)

        # We change the fontsize of minor ticks label
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.tick_params(axis='both', which='minor', labelsize=8)

        # Plot the usage of ECMs as vectors
        for i in range(len(Xusage)):
            color = cmap_curr[i, :]
            ax.quiver(Xusage[i], Yusage[i], Uusage[i], Vusage[i], pivot='tail', angles='xy', scale_units='xy',
                      linestyle='--', color=color, scale=1, width=quiverwidth)

        if objective_dict['obj_cons'] == 'obj':
            ax.set_title("Fraction needed \n per %.2f: " % scaling_factor + obj_string, size=8)
            # plt.title("Fraction needed per %.2f:\n" % scaling_factor + obj_string)
        else:
            ax.set_title("Fraction needed for \n demanded: " + obj_string, size=8)
            # plt.title("Fraction needed for demanded:\n" + obj_string)

        #ax.legend(loc='upper right')
        # set legend of axes invisible
        ax.get_legend().set_visible(False)
        # get legend handles and labels to create shared legend
        handles, labels = ax.get_legend_handles_labels()

    # Create common x and y axis labels
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)

    plt.xlabel('Needed fraction of constraint: ' + x_label)
    if not 'virtual' in cons_IDs:
        plt.ylabel(ylabel='Needed fraction of constraint: ' + y_label)

    # Create shared legend based on latest ax handles and labels
    fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.), loc='lower center',
               borderaxespad=0., ncol=2,
               frameon=False, prop={'size': 8})
    plt.subplots_adjust(bottom=0.1, hspace=0.6)

    if result_dir:
        fig.savefig(
            os.path.join(result_dir,
                         "cost_plot" + model_dict['model_name'] + "_" + cons_IDs[0] + "_" + cons_IDs[1] + ".png"),
            bbox_inches="tight")
    else:
        fig.savefig(
            os.path.join(os.getcwd(),
                         "cost_plot" + model_dict['model_name'] + "_" + cons_IDs[0] + "_" + cons_IDs[1] + ".png"),
            bbox_inches="tight")

def select_cons_IDs(infos_cons, cons_IDs):
    # If no cons_IDs are given, we take the first two in the list of constraints
    if not len(cons_IDs):
        cons_IDs = [infos_dict['ext_metab'] for infos_dict in infos_cons[:2]]
        print("No constraint-IDs were given. Therefore, the first two constraint shall be plotted, and the rest shall"
              "be ignored.")
    # If too many cons_IDs are given, we take only the first two of that list.
    elif len(cons_IDs) > 2:
        cons_IDs = cons_IDs[:2]
        print(
            "Too many constraint-IDs were given. Therefore, only the first two constraint-IDs shall "
            "be taken into account.")
    # If only one cons_ID is given. We will add a virtual constraint with zero costs.
    elif len(cons_IDs) == 1:
        cons_IDs.append('virtual')

    return cons_IDs


def get_cost_df(model_dict, infos_cons, cons_IDs, obj_index, obj_dict, scaling_factor=1):
    """
    This function takes in information about the ECMs, constraints and objectives, and produces a cost table for only
    the two constraintst that are to be plotted.
    If we consider a real objective, we rescale everything such that the costs are "the necessary fraction of the
    constraint to produce *scaling-factor* of objective flux.
    If we consider a secondary objective (a lower bound that has to be met) we rescale such that the costs are "the
    necessary fraction of the constraint to produce the needed amount of the secondary objective.
    :param model_dict: dictionary
            Information about the model for which the ECMs were calculated
    :param infos_cons: lists of dictionaries
            With information about the constraints
    :param cons_IDs: list of two strings
            Species-IDs of metabolites corresponding to the constraints for which the costs should be plotted
    :param obj_index: int
            The row at which the objective flux is found in the ecms_matrix
    :param obj_dict: python dictionary
            Dictionary with information about the objective
    :param scaling_factor:
            Scaling factor with which the costs will be scaled
    :return cost_df: Pandas dataframe
            Dataframe with for each ECM (columns) the costs with respect to the chosen constraints
    :return cons_indices: list of ints
            Indices at which, in the relevant_infos table, the chosen constraints can be found.
    """
    # Load relevant information
    original_table_cons_df = model_dict['table_cons_df']
    table_obj_df = model_dict['table_obj_df']
    # Make copy, because we do not want to affect the original dataframe when we are going to rescale
    table_cons_df = original_table_cons_df.copy()
    table_cons_df = table_cons_df.astype(float)

    # Loop over all constraints and drop it from the cost table if we do not want to plot this constraint now
    for cons_dict in infos_cons:
        cons_metab = cons_dict['ext_metab']
        if cons_metab not in cons_IDs:
            table_cons_df.drop(cons_metab, axis=1)

    # If we only have one real contraint, then the additional virtual constraint will have zero costs
    if 'virtual' in cons_IDs:
        table_cons_df['virtual'] = 0

    # Find indices at which, in the table_cons_df, the chosen constraints can be found.
    cons_indices = [index for index, column_label in enumerate(table_cons_df.columns) if
                    column_label in cons_IDs]
    # After this we have a dataframe with on each row an ECM, and as columns: constraint 1, constraint 2,
    # original_ECM_number, activity (optional) (but not necessarily in that order)

    # We will now normalize with respect to the objective
    # Find the column with the objective flux values for each ECM
    norm_column = table_obj_df.values[:, obj_index]
    # Drop zeros in normalization vector
    nonzeros = np.where(norm_column != 0)[0]
    norm_column = norm_column[nonzeros]
    # Determine norm_factor. This factor is used to get the costs per unit biomass production for the objective, and
    # the costs for satisfying the constraint for the lower bounds (that act as a kind of objective)
    if np.isnan(obj_dict['flux_val']):
        norm_factor = 1 * scaling_factor
    else:
        norm_factor = obj_dict['flux_val']
        # The abs in the following is necessary because the exchange reactions that are constrained are in fact removed
        # from the model. Therefore the metabolite seems to be taken up, while in fact the production of this metabolite
        # by the exchange reaction is constrained. I don't know how to fix this at the moment
        # TODO: Find something for this.
        norm_factor = norm_factor / abs(obj_dict['ext_metab_stoich'])

    # We want to find the costs as a fraction of the permitted flux through the constraint. So we have to divide the
    # costs by the bound for that constraint
    bound_values = [cons_dict['flux_val'] * abs(cons_dict['ext_metab_stoich']) for cons_dict in infos_cons if
                    cons_dict['ext_metab'] in cons_IDs]
    # The following is a (n_ECMs x 2) matrix with the appropriate factors for normalization
    norm_matrix = (np.transpose([norm_column]) * bound_values) / norm_factor

    # The activities should be scaled the other way as the costs, because the constraints should still be met with
    # equality
    if 'active' in table_cons_df.columns:
        norm_column_activities = norm_column / norm_factor

    # Make new dataframe in which zero rows are left out and is normalized
    cost_df = table_cons_df.copy()
    cost_df = cost_df.iloc[nonzeros, :]
    cost_array = cost_df.values[:, cons_indices]
    norm_cost_array = np.divide(cost_array, norm_matrix)
    cost_df.iloc[:, cons_indices] = norm_cost_array
    if 'active' in table_cons_df.columns:
        cost_df.loc[:, 'active'] = cost_df.loc[:, 'active'].multiply(norm_column_activities)

    return cost_df, cons_indices


def plot_different_approaches_one_figure(list_model_dicts, infos_obj, infos_cons, cons_IDs=[], obj_val=1,
                                         result_dir=None):
    """
    Plots the results of the different models (with different (sub)networks) in one plot.
    :param result_dir: directory path
            Directory for storing the figures
    :param list_model_dicts: list
            List of dictionaries with all information about one of the models that we are considering
    :param infos_obj: list
            List of dictionaries concerning the different objectives
    :param infos_cons: list
            List of dictionaries concerning the different constraints
    :param cons_IDs: list of cons_IDs
            List of strings reflecting the different constraints that should be plotted
    """
    # TODO: option to choose which constraints you want to show together, e.g. glc in combi with all other constraints

    # Some important numbers
    n_approaches = len(list_model_dicts)
    n_objectives = len(infos_obj)
    cmap = cm.get_cmap('tab10', n_approaches)
    scatter = np.linspace(10, 100, n_approaches)

    # We will use the scaling factor to get a nice overview of the plots. For example, when the obj_val is really large,
    # we will show the costs to produce more than one unit objective.
    scale_obj_val = np.floor(np.log2(obj_val))
    scaling_factor = 2 ** scale_obj_val

    # Make sure that we always have two constraints for the plot. 3_D plots are confusing
    # TODO: make a loop such that we can also plot all combinations of constraints ?
    cons_IDs = select_cons_IDs(infos_cons, cons_IDs)

    # Create one figure with n subplots, where n is the number of 'objectives'. A lower bound that should be met is also
    # considered an objective here.
    if n_objectives > 3:
        # Todo: change figsize dependent on number of (secondary) objectives
        fig, axes = plt.subplots(np.int(np.ceil(n_objectives/3)), 3, figsize=(7, 6))
    else:
        fig, axes = plt.subplots(1, n_objectives, figsize=(7, 4))
    if n_objectives == 1:
        axes = [axes]  # Necessary evil

    # We loop over the subplots, i.e., the number of objectives
    for index, objective_dict in enumerate(infos_obj):
        # Select the first subplot and the objective name
        if n_objectives > 3:
            ax = axes[np.int(index/3)][index%3]
        else:
            ax = axes[index]
        plt.sca(ax)
        obj_string = objective_dict['ext_metab']
        # In the x_ulims and y_ulims we will store a number of candidates for the upper bound of the plot. We will
        # eventually take the maximum of these candidates
        x_ulims = []
        y_ulims = []

        # Then loop over the different approaches that were taken, i.e., the different (sub)networks for which
        # cost-calculation was done
        for ind_approach in range(n_approaches - 1, -1, -1):
            # Select the dictionary with all the information
            model_dict = list_model_dicts[ind_approach]
            print()

            # filtered_cons_IDs = [rid for rid in cons_IDs if rid['rid'] in model_dict['model'].getReactionIDs()]
            filtered_cons_IDs = cons_IDs

            # Make cost table with only the constraints that are to be shown
            cost_df, cons_indices = get_cost_df(model_dict, infos_cons, filtered_cons_IDs, index, objective_dict,
                                                # was infos_cons, cons_IDs
                                                scaling_factor=scaling_factor)

            x_label = filtered_cons_IDs[0]  # cons_IDs[0]
            y_label = filtered_cons_IDs[1]  # cons_IDs[1]

            x_llim = y_llim = 0
            # Determine the bounds for the axes along which we want to plot the costs. It is a bit complicated, but
            # works
            # print(cons_indices)
            # print(cost_df)
            # print(max(cost_df.values[:, 0]))
            if objective_dict['obj_cons'] == 'obj':
                x_ulims.append(max(1.1, min(5 * scaling_factor / obj_val, max(cost_df.values[:, 0]) * 1.1)))
                y_ulims.append(max(1.1, min(5 * scaling_factor / obj_val, max(cost_df.values[:, 1]) * 1.1)))
            else:
                x_ulims.append(1.1)
                y_ulims.append(1.1)

            if 'virtual' in cons_IDs:
                y_ulims = [0.3]
                y_llim = -0.3

            # We want a different colour for the different approaches for which the results are plotted
            colors = np.tile(cmap.colors[ind_approach, :], (cost_df.shape[0], 1))

            # Make a scatter plot
            cost_df.plot.scatter(x=x_label, y=y_label, color=colors, ax=ax,
                                 s=scatter[ind_approach], label=list_model_dicts[ind_approach]['model_name'])

        # Set dimensions of the plot and configure axes
        # Also we make a light grey "constraint-box". The allowed region for solutions should be in this box.
        ax.set(adjustable='box', aspect='equal')
        y_ulim = max(y_ulims)
        x_ulim = max(x_ulims)
        plt.plot([1, 1], [y_llim, y_ulim], '--', color='grey', linewidth=2.0)
        ax.set_xlim([x_llim, x_ulim])
        ax.set_ylim([y_llim, y_ulim])
        plt.xlabel('') #'Needed fraction of constraint: ' + x_label)
        if 'virtual' in cons_IDs:
            ax.axes.get_yaxis().set_visible(False)
            ax.axhline(y=0, color='grey', linewidth=2.0)
        else:
            plt.ylabel('') #'Needed fraction of constraint: ' + y_label)
            plt.plot([0, 1], [0, 1], '-.', color='grey')
            ax.axes.get_yaxis().set_visible(True)
            plt.plot([0, x_ulim], [1, 1], '--', color='grey', linewidth=2.0)

        # We change the fontsize of minor ticks label
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.tick_params(axis='both', which='minor', labelsize=8)

        if objective_dict['obj_cons'] == 'obj':
            plt.title("Fraction needed \n per %.2f: " % scaling_factor + obj_string, size=8)
        else:
            plt.title("Fraction needed for \n demanded: " + obj_string, size=8)

        #plt.legend(loc='upper right')
        # set legend of axes invisible
        ax.get_legend().set_visible(False)
        # get legend handles and labels to create shared legend
        handles, labels = ax.get_legend_handles_labels()

    # Create shared x and y axis labels
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)

    plt.xlabel('Needed fraction of constraint: ' + x_label)
    if not 'virtual' in cons_IDs:
        plt.ylabel(ylabel='Needed fraction of constraint: ' + y_label)

    # Create shared legend based on latest ax handles and labels
    fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.), loc='lower center', borderaxespad=0.,
               ncol=np.int(np.ceil(n_approaches/2)) , frameon=False, prop={'size': 8})
    plt.subplots_adjust(bottom=0.15, hspace=0.70) #, wspace=0.3)

    if result_dir:
        # plt.tight_layout()
        plt.savefig(
            os.path.join(result_dir,
                         "various_approaches_constraints_" + cons_IDs[0] + "_" + cons_IDs[1] + ".png"))  # , dpi=600)
    else:
        # plt.tight_layout()
        plt.savefig(os.path.join(os.getcwd(), "various_approaches.png"))  # , dpi=600)


def find_hide_indices(model_path, to_be_tagged=[], focus_metabs=[], use_external_compartment=None):
    """
    :param model_path: string
    :param to_be_tagged: list of reaction-IDs
            List of reaction-IDs that should be tagged with a virtual metabolite
    :param focus_metabs:
            List of metabolites that cannot be ignored
    :return:
    """
    hide_indices = []
    if len(focus_metabs):
        # ECMtool works with a network class object. To find the indices of the metabolites that can be ignored, we need
        # this object already
        network = get_network_class(model_path, reactions_to_tag=to_be_tagged,
                                    use_external_compartment=use_external_compartment)
        for metab_index, metab in enumerate(network.metabolites):
            if metab.is_external:  # Internal metabolites are ignored by definition of conversion modes
                if metab.id not in focus_metabs:
                    hide_indices.append(metab_index)

    return hide_indices


def get_ecm_df(result_dir, model_path, infos_obj, to_be_tagged=[], print_results=False, hide_indices=[],
               get_EFMs=True, use_external_compartment=None, only_rays=False):
    """
    Gets cost table for each objective function.
    :param result_dir: string
    :param model_path: string
            Path to the xml-model for which ecms should be calculated
    :param infos_obj: list of dictionaries
            With information about the objectives
    :param to_be_tagged: list
            List of reactions that should be tagged with a virtual metabolite
    :param print_results: Boolean
            By default printing is off, because it can take a while with many ECMs
    :param hide_indices: list of ints
            By default, no metabolites are ignored. But we can do it
    :param external_compartment: string
            By default 'e', otherwise as given in the arguments with a string.
    :param get_EFMs: Boolean
            Calculate EFMs or not
    :return ecms_df: Pandas Dataframe
            Dataframe with as columns the ECMs and as rows the metabolite changes
    :return full_network_ecm: Network
            The network class for which the ecms were calculated
    :return efms_df: Pandas Dataframe
            Dataframe with as columns the different EFMs and as rows the fluxes
    :return ecms_from_efms: Pandas Dataframe
            Dataframe with all ECMs calculated from the EFMs by multiplying with stoichiometry
    """
    ecms_matrix, full_network_ecm = calc_ECMs(model_path, reactions_to_tag=to_be_tagged,
                                              print_results=print_results, hide_metabs=hide_indices,
                                              use_external_compartment=use_external_compartment, only_rays=only_rays)
    # print(ecms_matrix)

    # Normalize the ECMs. We first normalize all ECMs that produce the objective to produce one objective, after that
    # we normalize to different metabolite productions, such that the ECMs are maximally comparable
    ecms_matrix = normalize_ECMS_objective_first(ecms_matrix, full_network_ecm, infos_obj)

    # Find all metabolite_ids in the order used by ECMtool
    metab_ids = [metab.id for metab in full_network_ecm.metabolites]
    # Create dataframe with the ecms as columns and the metabolites as index
    ecms_df = pd.DataFrame(ecms_matrix, index=metab_ids)

    # Calculate EFMs if this is asked for. Infeasible for larger models.
    if get_EFMs:
        print("Calculating EFMs")
        full_network = copy.deepcopy(full_network_ecm)
        full_network.split_reversible()
        efms_df = calc_EFMs(full_network, result_dir)
        # Convert the EFMs to ECMs
        ecms_from_efms = convert_EFMs_to_ECMs(efms_df, full_network, infos_obj)
    else:
        ecms_from_efms = None
        efms_df = None

    return ecms_df, full_network_ecm, efms_df, ecms_from_efms


def check_bijection_Erik(ecms_first, ecms_second, network=None, normalization_order=[]):
    """
    :param ecms_first: np.array
            Matrix with ecms as columns, metabolites as rows
    :param ecms_second: np.array
            Matrix with ecms as columns, metabolites as rows
    :return bijection_YN: Boolean
    :return ecms_second_min_ecms_first: np.array
            Matrix with as columns the ECMs that were in the second but not in the first set
    :return ecms_first_min_ecms_second: np.array
            Matrix with as columns the ECMs that were in the first but not in the second set
    :param network: ecmtool.network class
            ecmtool.network object as comes from ECMtool and which is used for calculating ecms_matrix
    """
    # We first remove duplicates from both
    n_ecms_first_non_unique = ecms_first.shape[1]
    n_ecms_second_non_unique = ecms_second.shape[1]
    ecms_first = np.transpose(unique(np.transpose(ecms_first)))
    ecms_second = np.transpose(unique(np.transpose(ecms_second)))
    n_ecms_first = ecms_first.shape[1]
    n_ecms_second = ecms_second.shape[1]

    if n_ecms_first_non_unique - n_ecms_first > 0:
        print("Watch out. The first set of ECMs has duplicates")
    if n_ecms_second_non_unique - n_ecms_second > 0:
        print("Watch out. The second set of ECMs has duplicates")

    # Normalize both sets of ECMs
    if not len(normalization_order):
        normalization_order = determine_normalization_order(ecms_first, network)

    ecms_first = normalize_ECMS(ecms_first, network, normalization_order=normalization_order)
    ecms_second = normalize_ECMS(ecms_second, network, normalization_order=normalization_order)

    found_match_ecms_first = [False] * n_ecms_first
    no_match_ecms_second = list(range(n_ecms_second))
    for ecm_first_ind in range(n_ecms_first):
        if ecm_first_ind % 100 == 0:
            print('%d/%d ECMs checked for matches' % (ecm_first_ind, n_ecms_first))
        ecm_first = ecms_first[:, ecm_first_ind]
        for index, ecm_second_ind in enumerate(no_match_ecms_second):
            ecm_second = ecms_second[:, ecm_second_ind]

            if max(ecm_first - ecm_second) <= 10 ** -6:
                found_match_ecms_first[ecm_first_ind] = True
                del no_match_ecms_second[index]
                break

    ecms_first_min_ecms_second_inds = np.where([not found for found in found_match_ecms_first])[0]
    ecms_second_min_ecms_first_inds = no_match_ecms_second

    ecms_first_min_ecms_second = ecms_first[:, ecms_first_min_ecms_second_inds]
    ecms_second_min_ecms_first = ecms_second[:, ecms_second_min_ecms_first_inds]

    if not (ecms_first_min_ecms_second.shape[1] > 0 or ecms_second_min_ecms_first.shape[0] > 0):
        bijection_YN = True
    else:
        bijection_YN = False

    return bijection_YN, ecms_first_min_ecms_second, ecms_second_min_ecms_first


def determine_normalization_order(ecms_matrix, network):
    """
    Determine order of metabolites to which we are going to normalize.
    :return normalization_order: list of metab-IDs
            List of ordered metabolite-IDs
    :param ecms_matrix:
            This array contains the ECMs as columns and the metabolites as rows
    :param network: ecmtool.network class
            ecmtool.network object as comes from ECMtool and which is used for calculating ecms_matrix
    """
    metabs = []
    metabs_usage = []
    for metab_ind, metab in enumerate(network.metabolites):
        metabs.append(metab.id)  # Store all other metabolite ids
        metabs_usage.append(np.count_nonzero(ecms_matrix[metab_ind, :]))  # And their number of occurrences in ECMs

    # Order the tertiary objectives for how often they are used
    normalization_order = [x for (y, x) in
                           sorted(zip(metabs_usage, metabs), key=lambda pair: pair[0],
                                  reverse=True)]

    return normalization_order


def get_relevant_parts_ECMs(ecms_df, infos_cons, infos_obj):
    """
    Gets two dataframes: one with all the fluxes through all objectives for each ECM, one with all the fluxes through
    all constraints for each ECM. Note that we do not normalize in this function. The fluxes are equal to the ones in
    ecms_df, we just select the relevant information
    :param ecms_df: Pandas.Dataframe
            A cost_df with ECMs as columns and as rows the metabolites
    :param infos_cons: list of dictionaries
            List of dictionaries with information on all active constraints
    :param infos_obj: list of dictionaries
            List of dictionaries with information on all active objectives (including secondary objectives such as
            lower bounds that have to be met.
    """
    # Make local copy, because I don't want to change the original ecms
    ecms_df_loc = ecms_df.copy()

    # Initialize some stuff
    n_ECMs = ecms_df_loc.shape[1]
    n_obj = len(infos_obj)
    n_cons = len(infos_cons)
    table_cons = np.ones((n_cons, n_ECMs))
    table_obj = np.ones((n_obj, n_ECMs))
    names_obj = []
    names_cons = []
    counter_obj = 0
    counter_cons = 0

    # Loop over all objectives and constraints
    for obj_cons_dict in infos_obj + infos_cons:
        # Get Species-ID for metabolite that is coupled to the current objective/constraint
        sid_cons_obj = obj_cons_dict['ext_metab']

        # Get production/consumption flux for above metabolite
        obj_cons_flux = ecms_df_loc.loc[sid_cons_obj, :].values

        # Add the flux information to a table
        if obj_cons_dict['obj_cons'] == 'cons':
            table_cons[counter_cons, :] = obj_cons_flux
            names_cons.append(sid_cons_obj)
            counter_cons += 1
        else:
            table_obj[counter_obj, :] = obj_cons_flux
            names_obj.append(sid_cons_obj)
            counter_obj += 1

    table_obj_df = pd.DataFrame(table_obj.transpose(), columns=names_obj)
    table_cons_df = pd.DataFrame(table_cons.transpose(), columns=names_cons)

    # Add number-id to ECMs. Needed for when we start deleting inactive ECMs and so on.
    orig_ECM_number = list(range(n_ECMs))
    table_obj_df['orig_ECM_number'] = orig_ECM_number
    table_cons_df['orig_ECM_number'] = orig_ECM_number

    return table_cons_df, table_obj_df


def convert_flux_vector_to_conversion(flux_vector, cmod, network):
    """
    Gets activities of different EFMs according to FBA solution.
    :param flux_vector: ndarray
    :param cmod: cbmpy.CBmodel
    :param network: ecmtool.network class
    :return conv_vector: ndarray
    """

    # Delete exchange reactions from the flux vector, so that we can calculate the conversion
    pairs = cbm.CBTools.findDeadEndReactions(cmod)
    reac_indices = []
    external_metabolites, external_reactions = zip(*pairs) if len(pairs) else (
        zip(*cbm.CBTools.findDeadEndMetabolites(cmod))[0], [])
    for re_index, reac in enumerate(cmod.reactions):
        re_id = reac.id
        if re_id in external_reactions:
            reac_indices.append(re_index)

    # Find flux vector for only the internal fluxes
    flux_vector_int = flux_vector[np.setdiff1d(range(len(flux_vector)), reac_indices)]

    conv_vector = np.dot(network.N, flux_vector_int)

    return conv_vector


# def unique(matrix): ## eunice_edit
#    unique_set = {tuple(row) for row in matrix if np.count_nonzero(row) > 0}
#    return np.vstack(unique_set) if len(unique_set) else to_fractions(np.ndarray(shape=(0, matrix.shape[1])))


def convert_EFMs_to_ECMs(efms_df, network, infos_obj, verbose=True):
    """
    Gets activities of different EFMs according to FBA solution.
    :param infos_obj: list of dictionaries
            With information needed for normalization with respect to objectives
    :param efms_df: pandas Dataframe
            EFMs as columns, reactions as rows
    :param network: ecmtool.network class
    :param verbose: Boolean
    :return ecms_df: pandas dataframe
    """

    n_EFMs = efms_df.shape[1]
    if verbose:
        print("Converting EFMs to ECMs.")

    # Select only the external metabolites, the internal metabolite changes are zero anyhow in the EFMs
    ext_indices = [index for index, metabolite in enumerate(network.metabolites) if metabolite.is_external]
    ext_stoichiometry = network.N[ext_indices, :]

    ecms_array = np.zeros((network.N.shape[0], n_EFMs))
    # Multiply EFMs with the stoichiomery per 100 EFMs. Otherwise, it can be too big a task.
    for i in range(int(n_EFMs / 100)):
        if verbose:
            print('%d/%d have been converted' % (100 * i, n_EFMs))
        flux_vectors = efms_df.values[:, 100 * i:100 * (i + 1)]
        conv_vector = np.dot(ext_stoichiometry, flux_vectors)
        ecms_array[ext_indices, 100 * i:100 * (i + 1)] = conv_vector

    # And do the remaining EFMs. Of course this could have been done more elegantly
    i = int(n_EFMs / 100)
    if verbose:
        print('%d/%d have been converted' % (100 * i, n_EFMs))
    flux_vectors = efms_df.values[:, 100 * i:]
    conv_vector = np.dot(ext_stoichiometry, flux_vectors)
    ecms_array[ext_indices, 100 * i:] = conv_vector
    if verbose:
        print('%d/%d have been converted' % (n_EFMs, n_EFMs))

    if verbose:
        print("Normalizing converted EFMs")
    ecms_array = normalize_ECMS_objective_first(ecms_array, network, infos_obj)

    # Round off to some decimal places. It seems that EFMtool introduces many round-off errors.
    ext_ecms_array = ecms_array[ext_indices, :]
    # print(ext_ecms_array)
    # print(np.size(ecms_array))
    # print(type(ext_ecms_array))
    # ecms_array[ext_indices, :] = np.around(ext_ecms_array.astype(np.double), decimals=2)
    ext_ecms_array_rounded = np.around(ext_ecms_array.astype(np.double), decimals=2)

    # Keep only the unique ECMs
    # ecms_array = np.transpose(unique(np.transpose(ecms_array)))
    dummy, unique_inds = np.unique(ext_ecms_array_rounded, axis=1, return_index=True)
    ext_ecms_array = ext_ecms_array[:, unique_inds]

    metab_ids = [metab.id for metab in network.metabolites if metab.is_external]
    ecms_df = pd.DataFrame(ext_ecms_array, index=metab_ids)
    if verbose:
        print('Found ' + str(ext_ecms_array.shape[1]) + ' unique conversions generated by the ' + str(
            n_EFMs) + ' EFMs')
    return ecms_df


def get_FBA_result(cmod, network):
    """
    Uses the cbmpy model and the ecmtool.network class to find for each reaction in network the corresponding flux in
    the FBA solution.
    Note: Some reactions in the cbmpy model are not in the network class. Mostly these are exchange reactions. These are
    deleted in ecmtool, because otherwise no conversions can be calculated (everything will be in steady state)
    :param cmod: cbmpy model object
            FBA should already have been done
    :param network: ecmtool.network object
    :return: FBA_vector: np.array
            Vector with for each reaction in network the FBA flux through it
    """
    # Get list of reactions in the network object, and initialize a zero vector
    reac_ids = [reac.id for reac in network.reactions]
    FBA_vector = np.zeros((len(reac_ids), 1))

    # Loop over all reactions in network, find their cbmpy equivalent and store their FBA-value
    for index, reac_id in enumerate(reac_ids):
        if reac_id in cmod.getReactionIds():
            cmod_rea = cmod.getReaction(reac_id)
            FBA_val = cmod_rea.getValue()
        else:  # This reaction does not exist in cmod
            FBA_val = 0

        FBA_vector[index] = FBA_val

    return FBA_vector


def get_activity_ECMs(ecms_df, cmod, network, round_off_to_zero=True, hide_indices=[]):
    """
    Gets activities of different EFMs according to FBA solution.
    :param hide_indices: list of ints
            List of metabolite indices that were ignored in emc-calculation and should thus also be ignored in fitting
            the FBA solution
    :param round_off_to_zero: Boolean
            Boolean setting if small activities of fluxes can be rounded off
    :param network: ecmtool.network class
            The network object for which ecms_df was calculated
    :param ecms_df: Pandas.Dataframe
            A cost_df with ECMs as columns and as rows the metabolites
    :param cmod: cbmpy.CBModel
            An FBA should already have been performed
    :return ecm_usage: np.array
            Vector of length the number of ECMs with the activities that the ECMs should have to fit the FBA solution
    """
    # TODO: Somehow this function doesn't really work when there is a lowerbound constraint (kind of additional objective)
    # Do you mean a lowerbound constraint >0 ? or an upperbound constraint <0. This resulted in resp negative and positive flux with
    # a resp positive and negative nzrc.
    # If the sign of flux value and nzrc is equal it becomes a constraint, if different it becomes a (secondary) objective.
    # In the latter case the flux should become smaller (i.e. closer to zero) to improve the objective.
    # In the first case (flux value and nzrc with equal sign) the flux should become bigger (more negative or more positive).

    # We still need to find out why this -the problem with additional objectives- happens

    LP_TOLERANCE = 10 ** -7 #9 #10
    ZERO_TOLERANCE = 10 ** -7 #12 was 8

    # Get flux vector solution from FBA
    flux_solution = get_FBA_result(cmod, network)

    # Convert it to a conversion
    conv_solution = np.dot(network.N, flux_solution)

    ecms_array = ecms_df.values  # The rows are the metabolites, columns are conversions

    # If the ECMs are calculated while ignoring some of the external metabolites, then these should also be ignored
    # when we try to fit the FBA solution with the calculated ECMs
    to_be_deleted = []
    if len(hide_indices):
        to_be_deleted = to_be_deleted + hide_indices

    # Check if internal metabolites are indeed (close to) zero. If so, throw them away before doing the LP, because they
    # will only lead to redundancy in the LP
    internal_inds = [ind for ind, metab in enumerate(network.metabolites) if not metab.is_external]
    internal_not_hidden = [ind for ind in internal_inds if ind not in hide_indices]
    if np.max(np.abs(conv_solution[internal_not_hidden])) > ZERO_TOLERANCE:
        print('Warning: some internal metabolites seem not to be in steady state in the FBA solution')
    else:
        to_be_deleted = to_be_deleted + internal_inds

    to_be_deleted = np.unique(to_be_deleted)

    # Constraint matrix is A_eq: rows are the external (not-ignored) metabolites. For each metabolite the weighted
    # sum of ECMs should be equal to the conversion derived from the FBA-solution (captured by b_eq)
    A_eq = np.delete(ecms_array, to_be_deleted, axis=0)
    b_eq = np.delete(conv_solution, to_be_deleted)
    # As an objective, we try to minimize the sum of weights. This shouldn't really matter for the result.
    c = np.ones((A_eq.shape[1], 1))

    result = linprog(c, A_eq=A_eq, b_eq=b_eq,
                     method='simplex', options={'tol': LP_TOLERANCE})

    if not result.success:
        print('Couldn\'t fit the FBA solution, so no activities are calculated')
    else:
        # TODO: Remove the following print-statement or add annotation
        print(result.fun)

    ecm_usage = result.x.transpose()

    # The LP-solver has some tolerance so that some weights will be non-zero. We here round these off to zero.
    if round_off_to_zero:
        lower_bound = ZERO_TOLERANCE * ecm_usage.max()
        nonzeros_before = np.count_nonzero(ecm_usage)
        ecm_usage[np.where(ecm_usage < lower_bound)[0]] = 0
        nonzeros_after = np.count_nonzero(ecm_usage)
        n_rounded_off = nonzeros_before - nonzeros_after
        if n_rounded_off:
            print(
                'Warning: The activities of %d ECMs are being round off to zero.\n '
                'Select round_off_to_zero = False if you don\'t want this.' % n_rounded_off)

    # Check if this is equal to what we expect
    predicted_conv = np.dot(ecms_array, ecm_usage)
    difference = abs(np.delete(predicted_conv - np.transpose(conv_solution), hide_indices))
    max_diff = difference.max()
    print('Difference in the supremum-norm between estimated ECMs and FBA solution is: %f' % max_diff)

    return ecm_usage


# TODO: This function is copied from ecmtool. Check if this is okay.
def print_ECMs(cone, debug_tags, network, orig_N, add_objective_metabolite, check_feasibility):
    print("Printing ECM results")
    for index, ecm in enumerate(cone):
        # Normalise by objective metabolite, if applicable
        objective_index = -1 - len(debug_tags)
        objective = ecm[objective_index]
        if add_objective_metabolite and objective > 0:
            ecm /= objective

        metabolite_ids = [met.id for met in
                          network.metabolites] if not network.compressed else network.uncompressed_metabolite_ids

        print('\nECM #%d:' % index)
        for metabolite_index, stoichiometry_val in enumerate(ecm):
            if stoichiometry_val != 0.0:
                print('%s\t\t->\t%.4f' % (metabolite_ids[metabolite_index], stoichiometry_val))

        # if check_feasibility:
        #     allowed_error = 10 ** -6
        #     solution = linprog(c=[1] * orig_N.shape[1], A_eq=orig_N, b_eq=cone[index, :],
        #                        bounds=[(-1000, 1000)] * orig_N.shape[1], options={'tol': allowed_error})
        #     print('ECM satisfies stoichiometry' if solution.status == 0 else 'ECM does not satisfy stoichiometry')


def check_bijection_csvs(ecms_first_df, ecms_second_df):
    """
    :param ecms_first_df: DataFrame
            Matrix with ecms as rows, metabolites as cols
            colname should give metab_id
    :param ecms_second_df: DataFrame
            Matrix with ecms as rows, metabolites as cols
            colname should give metab_id
    :return bijection_YN: Boolean
    :return ecms_second_min_ecms_first: np.array
            Matrix with as columns the ECMs that were in the second but not in the first set
    :return ecms_first_min_ecms_second: np.array
            Matrix with as columns the ECMs that were in the first but not in the second set
    """
    # We first remove duplicates from both
    metab_ids_first = list(ecms_first_df.columns)
    metab_ids_second = list(ecms_second_df.columns)
    ecms_first = np.transpose(ecms_first_df.values)
    ecms_second = np.transpose(ecms_second_df.values)
    n_ecms_first_non_unique = ecms_first.shape[1]
    n_ecms_second_non_unique = ecms_second.shape[1]
    ecms_first = np.transpose(unique(np.transpose(ecms_first)))
    ecms_second = np.transpose(unique(np.transpose(ecms_second)))
    n_ecms_first = ecms_first.shape[1]
    n_ecms_second = ecms_second.shape[1]

    # Find matching of metab_ids
    matching_inds = np.zeros(len(metab_ids_first))
    for id_ind, id in enumerate(metab_ids_first):
        matching_inds[id_ind] = [id_ind_sec for id_ind_sec, id_sec in enumerate(metab_ids_second) if id_sec == id][0]

    # Make sure that second ecms metabolites are in the same order
    matching_inds = matching_inds.astype(int)
    ecms_second = ecms_second[matching_inds, :]

    if n_ecms_first_non_unique - n_ecms_first > 0:
        print("Watch out. The first set of ECMs has duplicates")
    if n_ecms_second_non_unique - n_ecms_second > 0:
        print("Watch out. The second set of ECMs has duplicates")

    # Normalize both sets of ECMs
    sum_columns_first = np.sum(np.abs(ecms_first), axis=0)
    sum_columns_first = sum_columns_first[np.newaxis, :]
    ecms_first = ecms_first / np.repeat(sum_columns_first, ecms_first.shape[0], axis=0)

    sum_columns_second = np.sum(np.abs(ecms_second), axis=0)
    sum_columns_second = sum_columns_second[np.newaxis, :]
    ecms_second = ecms_second / np.repeat(sum_columns_second, ecms_second.shape[0], axis=0)

    found_match_ecms_first = [False] * n_ecms_first
    no_match_ecms_second = list(range(n_ecms_second))
    for ecm_first_ind in range(n_ecms_first):
        if ecm_first_ind % 100 == 0:
            print('%d/%d ECMs checked for matches' % (ecm_first_ind, n_ecms_first))
        ecm_first = ecms_first[:, ecm_first_ind]
        for index, ecm_second_ind in enumerate(no_match_ecms_second):
            ecm_second = ecms_second[:, ecm_second_ind]

            if max(ecm_first - ecm_second) <= 10 ** -3:
                found_match_ecms_first[ecm_first_ind] = True
                del no_match_ecms_second[index]
                break

    ecms_first_min_ecms_second_inds = np.where([not found for found in found_match_ecms_first])[0]
    ecms_second_min_ecms_first_inds = no_match_ecms_second

    ecms_first_min_ecms_second = ecms_first[:, ecms_first_min_ecms_second_inds]
    ecms_second_min_ecms_first = ecms_second[:, ecms_second_min_ecms_first_inds]

    if not (ecms_first_min_ecms_second.shape[1] > 0 or ecms_second_min_ecms_first.shape[1] > 0):
        bijection_YN = True
    else:
        bijection_YN = False

    return bijection_YN, ecms_first_min_ecms_second, ecms_second_min_ecms_first


def get_active_ecms(model_dict):
    # TODO: find active ecms and select only the relevant input output metabolites with usage values compared ...
    # to one unit objective

    # Todo: add (secondary) objectives to this list.
    idx = list(model_dict['table_cons_df']['active'][model_dict['table_cons_df']['active'] != 0.].index)

    #model_dict['ecms_df'][model_dict['ecms_df'][idx] != 0]

    a = model_dict['ecms_df'][idx] != 0

    relevant_part_active_ecms = model_dict['ecms_df'][idx].loc[a.sum(axis=1) > 0]
    #print(relevant_part_active_ecms)
    #print(relevant_part_active_ecms.astype(float))

    return relevant_part_active_ecms


def delete_bounds_from_model(cmod):
    cmod_wo_bounds = cmod.clone()
    all_bounds = cmod.getAllFluxBounds()
    max_bound = max(abs(np.array(list(all_bounds.values()))))
    for bound in all_bounds:
        bound_value = all_bounds[bound]
        if (bound_value != 0.0) & (abs(bound_value) < max_bound):
            if bound[-3:] == '_lb':
                if bound_value < 0.0:
                    new_bound_value = -max_bound
                else:
                    new_bound_value = 0.0
            elif bound[-3:] == '_ub':
                if bound_value < 0.0:
                    new_bound_value = 0.0
                else:
                    new_bound_value = max_bound
            else:
                raise NameError('This does not work, because bound-IDs are not correct.')

            all_bounds[bound] = new_bound_value

    cmod_wo_bounds.setFluxBoundsFromDict(all_bounds)
    return cmod_wo_bounds


def find_associated_efms(cmod, table_cons_df, ecms_df, infos_obj_cons, model_path, use_external_compartment=None,
                         only_relevant_ECMs=True):
    """
    :param cmod: cbmpy.CBModel
    :param table_cons_df: pandas.Dataframe
            Dataframe with all the fluxes through all constraints for each ECM, and a column 'active' that
            indicates the activity of each ECM in the FBA solution.
    :param ecms_df: Pandas.Dataframe
            A cost_df with ECMs as columns and as rows the metabolites
    :param infos_obj_cons: list of dictionaries
            List of dictionaries with information on all active constraints and on all active objectives
            (including secondary objectives such as lower bounds that have to be met.
    :param model_path
    :param use_external_compartment: string
            If given, defines the way an external compartment is notated in the model (cmod).
    :param only_relevant_ECMs: Boolean
            If true only the active ECMs are used to find corresponding EFMs and full ECMs
            If false all ECMs are used.
    :return relevant_efms_df: pandas.Dataframe
            Matrix with as columns the reaction IDs and as index ECM IDs
    :return full_relevant_ecms_df: pandas.Dataframe
            Matrix with as columns the external metabolite IDs and as index ECM IDs
    """
    ZERO_TOLERANCE = 10 ** -6 # 1e-10

    # First delete old bounds from the model
    network = extract_sbml_stoichiometry(model_path, determine_inputs_outputs=True,
                                         use_external_compartment=use_external_compartment)
    cmod_wo_bounds = delete_bounds_from_model(cmod)
    all_bounds = cmod_wo_bounds.getAllFluxBounds()
    if not only_relevant_ECMs: #use all ECMs
        relevant_ecms_df = ecms_df.iloc[np.where(np.count_nonzero(ecms_df.values, axis=1) > 0)[0], :]
    else: #use only the active ECMs
        active_ECMs = table_cons_df[table_cons_df['active'] != 0.]['orig_ECM_number'].values
        relevant_ecms_df = ecms_df.iloc[np.where(np.count_nonzero(ecms_df.values, axis=1) > 0)[0], active_ECMs]
    all_rids = cmod_wo_bounds.getReactionIds()
    relevant_efms_vals = np.zeros((relevant_ecms_df.shape[1], len(all_rids)))
    external_metabs = [ind for ind, metab in enumerate(network.metabolites) if metab.is_external]
    full_ecms_vals = np.zeros((relevant_ecms_df.shape[1], len(external_metabs)))

    # Then loop over the different ECMs of interest
    for ecm_ind in range(relevant_ecms_df.shape[1]):
        ecm = relevant_ecms_df.iloc[:, ecm_ind]
        for metab_key in ecm.keys():
            info_metab = [dictionary for dictionary in infos_obj_cons if dictionary['ext_metab'] == metab_key][0]
            rid = info_metab['rid']
            metab_prod = float(ecm[metab_key])
            flux_value = metab_prod / abs(info_metab['ext_metab_stoich'])

            # Then set bounds that are given by the ECM of interest
            all_bounds[rid+'_lb'] = flux_value
            all_bounds[rid+'_ub'] = flux_value

        cmod_wo_bounds.setFluxBoundsFromDict(all_bounds)
        cbm.doFBA(cmod_wo_bounds)
        cbm.CBCPLEX.cplx_MinimizeNumActiveFluxes(cmod_wo_bounds)

        corr_flux_vals = [cmod_wo_bounds.getReaction(rid).value for rid in all_rids]
        relevant_efms_vals[ecm_ind, :] = corr_flux_vals

        # TODO: Find conversion by multiplying with network.N
        # Get list of reactions in the network object, and initialize a zero vector
        network_rids = [reac.id for reac in network.reactions]
        corr_ecm_vals = [cmod_wo_bounds.getReaction(rid).value for rid in network_rids]
        corr_conversion = np.dot(network.N, corr_ecm_vals)

        # Check if internal metabolites are indeed (close to) zero. If so, throw them away before doing the LP, because they
        # will only lead to redundancy in the LP
        internal_inds = [ind for ind, metab in enumerate(network.metabolites) if not metab.is_external]
        if np.max(np.abs(corr_conversion[internal_inds])) > ZERO_TOLERANCE:
            print(np.max(np.abs(corr_conversion[internal_inds])))
            raise ValueError('Warning: some internal metabolites seem not to be in steady state in the FBA solution')

        clean_conversion = np.delete(corr_conversion, internal_inds)
        print(ecm_ind)
        full_ecms_vals[ecm_ind, :] = clean_conversion

    if not only_relevant_ECMs:
        relevant_efms_df = pd.DataFrame(relevant_efms_vals, columns=all_rids, index=ecms_df.index)
        external_metab_ids = [network.metabolites[ind].id for ind in external_metabs]
        full_relevant_ecms_df = pd.DataFrame(full_ecms_vals, columns=external_metab_ids, index=ecms_df.index)
    else: # give index to indicate the corresponding ECM number
        relevant_efms_df = pd.DataFrame(relevant_efms_vals, columns=all_rids, index=active_ECMs)
        external_metab_ids = [network.metabolites[ind].id for ind in external_metabs]
        full_relevant_ecms_df = pd.DataFrame(full_ecms_vals, columns=external_metab_ids, index=active_ECMs)
    return relevant_efms_df, full_relevant_ecms_df


def get_cons_ID_combinations(cons_IDs):
    """Create list of lists of two cons IDs that will be plotted in a  cost vector graphs."""
    cons_ID_combi_list = []
    for i in range(int(np.ceil(len(cons_IDs) / 2))):
        cons_ID_combi_list = cons_ID_combi_list + [cons_IDs[i * 2: (i + 1) * 2]]
    return cons_ID_combi_list


def get_costs_ECMs(model_dict, infos_cons, cons_IDs):
    """ find needed fraction of constraints per ECM """

    # Load relevant information
    original_table_cons_df = model_dict['table_cons_df']
    table_obj_df = model_dict['table_obj_df']
    # Make copy, because we do not want to affect the original dataframe when we are going to rescale
    table_cons_df = original_table_cons_df.copy()
    table_cons_df = table_cons_df.astype(float)

    # Loop over all constraints and drop it from the cost table if we do not want to plot this constraint now
    for cons_dict in infos_cons:
        cons_metab = cons_dict['ext_metab']
        if cons_metab not in cons_IDs:
            table_cons_df.drop(cons_metab, axis=1)

    # If we only have one real contraint, then the additional virtual constraint will have zero costs
    if 'virtual' in cons_IDs:
        table_cons_df['virtual'] = 0

    # Find indices at which, in the table_cons_df, the chosen constraints can be found.
    cons_indices = [index for index, column_label in enumerate(table_cons_df.columns) if
                    column_label in cons_IDs]
    # After this we have a dataframe with on each row an ECM, and as columns: constraint 1, constraint 2,
    # original_ECM_number, activity (optional) (but not necessarily in that order)

    # We want to find the costs as a fraction of the permitted flux through the constraint. So we have to divide the
    # costs by the bound for that constraint
    # cons_IDs = ['M_glc__D_e', 'M_o2_e']
    norm_list = []
    for cons_ID in cons_IDs:
        # find bound
        total = sum(model_dict['table_cons_df']['active'] * model_dict['table_cons_df'][cons_ID])
        # print(model_dict['table_cons_df']['active']*model_dict['table_cons_df'][cons_ID]/total)
        norm_list.append(total)
    norm_matrix = np.array(norm_list)

    # model_dict['table_cons_df']['active']*model_dict['table_cons_df'][cons_ID]/10.

    # Make new dataframe in which zero rows are left out and is normalized
    cost_df = table_cons_df.copy()
    #cost_df = cost_df.iloc[nonzeros, :]
    cost_array = cost_df.values[:, cons_indices]
    #  divide by total available flux
    norm_cost_array = np.divide(cost_array, norm_matrix)
    # multiply with activities
    #norm_cost_array = norm_cost_array*np.transpose(np.array([cost_df['active'].values, cost_df['active'].values]))
    cost_df.iloc[:, cons_indices] = norm_cost_array
    
    return cost_df, cons_indices


def plot_costs_ECMs(model_dict, infos_cons, cons_IDs=[], show_active=True, result_dir=None, vector_focus=True):
    """
    Plots the costs of all ECMs in cost_df in the cost vector figure, and shows ECM_usage in the FBA-solution by having
    vectors.
    Let them sum up to 1. So not (secundary) objective specific. But will always end up in left top corner.
    :param result_dir: directory path
            Directory for storing the figures
    :param model_dict: dictionary
            Dictionary with all information about one of the models that we are considering
    :param infos_cons: list
            List of dictionaries concerning the different constraints
    :param cons_IDs: list of cons_IDs
            List of strings reflecting the different constraints that should be plotted
    :param show_active: Boolean
            Show or not show ECM usage
    :param vector_focus: Boolean
            Focus on vectors or show all ECMs
    """
    # Todo: scale ECMs such that they are plotted in the visible area.
    # Todo: scale based on objective? or such that at least one of the constraint values is plotted between 0 and 1.1 or 1.5?
    # Todo: adjust this at least for the active ECMs
    # Todo: use plot_various_approaches ?

    # Some important numbers
    #n_objectives = len(infos_obj)
    scatter_active = 100
    alpha_active = 0.5
    scatter_normal = 15
    quiverwidth = 0.01

    # If the objective is really large, we will show the costs to produce more than one unit objective. We use the
    # scaling factor to calculate the new costs
    #scale_obj_val = np.floor(np.log2(obj_val))
    #scaling_factor = 2 ** scale_obj_val

    # Find the indices of the ECMs that are active in the FBA solution
    actECMs_inds = list(model_dict['table_cons_df'].loc[model_dict['table_cons_df']['active'] > 0]['orig_ECM_number'])

    # We use gray for normal non-active ECMs, and a different colour for each active ECM
    # Todo: increase list of colors,
    # e.g. https://stackoverflow.com/questions/4971269/how-to-pick-a-new-color-for-each-plotted-line-within-a-figure-in-matplotlib
    cmap = cm.get_cmap('tab10', len(actECMs_inds))


    # Make sure that we always have two constraints for the plot
    cons_IDs = select_cons_IDs(infos_cons, cons_IDs)

    # Create one figure with n subplots, where n is the number of 'objectives'. A lower bound that should be met is also
    # considered an objective here.
    #if n_objectives > 3:
    #    fig, axes = plt.subplots(np.int(np.ceil(n_objectives/3)), 3, figsize=(7,6))
    #else:
    #    fig, axes = plt.subplots(1, n_objectives, figsize=(7, 4))
    fig, axes = plt.subplots(1, 1, figsize=(5,4))
    axes = [axes]  # Necessary evil
    # We loop over the subplots, i.e., the number of objectives
    ax = axes[0]
    plt.sca(ax)
    #obj_string = objective_dict['ext_metab']

    # In the x_ulims and y_ulims we will store a number of candidates for the upper bound of the plot. We will
    # eventually take the maximum of these candidates
    x_ulims = []
    y_ulims = []
    x_llims = []
    y_llims = []

    # Make cost table with only the constraints that are to be shown
    cost_df, cons_indices = get_costs_ECMs(model_dict, infos_cons, cons_IDs)
    print(cost_df)

    if show_active and 'active' not in cost_df.columns:
        print("Can't show active EFMs if this information is not provided.")
        return

    x_label = cons_IDs[0]
    y_label = cons_IDs[1]

    # Determine the bounds for the axes along which we want to plot the costs. It is a bit complicated, but
    # works

    x_ulims.append(1.1)
    y_ulims.append(1.1)
    x_llims.append(-0.05)
    if not vector_focus:
        x_ulims.append(max(cost_df[cons_IDs[0]]) * 1.1)
        x_llims.append(min(0, min(cost_df[cons_IDs[0]])))
    y_llims.append(-0.05)
    if not vector_focus:
        y_ulims.append(max(cost_df[cons_IDs[1]]) * 1.1)
        y_llims.append(min(0, min(cost_df[cons_IDs[1]])))

    if 'virtual' in cons_IDs:
        y_ulims = [0.3]
        y_llims = [-0.3]

    # Plot cost dots of normal ECMs
    scaled_cost_df = cost_df.copy()
    cost_array = cost_df.values[:, cons_indices]
    #  divide by 2
    norm_matrix = np.ones(np.shape(cost_array))*2
    norm_cost_array = np.divide(cost_array, norm_matrix)
    # multiply with activities
    # norm_cost_array = norm_cost_array*np.transpose(np.array([cost_df['active'].values, cost_df['active'].values]))
    scaled_cost_df.iloc[:, cons_indices] = norm_cost_array
    scaled_cost_df['active'] = np.multiply(scaled_cost_df['active'], 2)
    scaled_cost_df.plot.scatter(x=x_label, y=y_label, color='grey', ax=ax,
                         s=scatter_normal, label=model_dict['model_name']) #, legend=False)

    # Plot cost dots of active ECMs
    actECMs = scaled_cost_df.loc[scaled_cost_df['orig_ECM_number'].isin(actECMs_inds)]
    # The number of active ECMs for this objective can be lower than the number of active ECMs in total, because
    # we cannot plot the costs for an ECM that does not contribute anything to this objective.
    n_actECMs_curr_obj = len(actECMs)

    # Here we will store the x,y-length of the cost vectors if multiplied by its activity
    ecm_usage = np.zeros((n_actECMs_curr_obj, 2))
    # x-length
    ecm_usage[:, 0] = actECMs.values[:, cons_indices[0]] * actECMs['active']
    # y-length
    ecm_usage[:, 1] = actECMs.values[:, cons_indices[1]] * actECMs['active']

    # Construct two columns of starting positions for the vectors for the vector addition with which we show the
    # usage of the ECMs. The end point of the first vector is the starting point of the second vector, etc.
    quiver_start = ecm_usage.cumsum(axis=0)
    # Append starting point (0,0) of the first vector
    quiver_start = np.append(np.zeros((1, 2)), quiver_start, axis=0)
    # Delete last row because this is actually the end point of the last vector
    quiver_start = quiver_start[:-1, :]
    # Add the x,y-lengths of the vectors as two additional columns
    ecm_usage = np.append(quiver_start, ecm_usage, axis=1)

    # The following thus are the x,y-begin positions of the vectors, and the x,y-length
    Xusage, Yusage, Uusage, Vusage = zip(*ecm_usage)

    # Get the indices of the ECMs that contribute to this objective, needed for getting the right colour
    curr_act_inds = [counter for counter, ind in enumerate(actECMs_inds) if ind in actECMs['orig_ECM_number']]
    cmap_curr = cmap.colors[curr_act_inds, :]

    # Plot the costs per unit objective of the active ECMs
    actECMs.plot(kind='scatter', x=x_label, y=y_label, color=cmap_curr, ax=ax,
                 s=scatter_active, alpha=alpha_active, label='active ECM', legend=True) #, label=[str(i) for i in actECMs_inds])

    # Plot the usage of ECMs as vectors
    # for i in range(len(Xusage)):
    #    color = cmap_curr[i, :]
    #    ax.quiver(Xusage[i], Yusage[i], Uusage[i], Vusage[i], pivot='tail', angles='xy', scale_units='xy',
    #              linestyle='--', color=color, scale=1, width=quiverwidth)

    # Set dimensions of the plot and configure axes
    # Also we make a light grey "constraint-box". The allowed region for solutions should be in this box.
    if vector_focus:
        ax.set(adjustable='box', aspect='equal')
    y_ulim = max(y_ulims)
    x_ulim = max(x_ulims)
    x_llim = min(x_llims)
    y_llim = min(y_llims)
    ax.plot([1, 1], [y_llim, y_ulim], '--', color='grey', linewidth=2.0)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_ylim([y_llim, y_ulim])
    # plt.xlabel('Needed fraction of constraint: ' + x_label)
    ax.set(xlabel='')
    if 'virtual' in cons_IDs:
        ax.axes.get_yaxis().set_visible(False)
        ax.axhline(y=0, color='grey', linewidth=2.0)
    else:
        ax.set(ylabel='')
        # plt.ylabel('Needed fraction of constraint: ' + y_label)
        ax.plot([0, 1], [0, 1], '-.', color='grey')
        ax.axes.get_yaxis().set_visible(True)
        ax.plot([x_llim, x_ulim], [1, 1], '--', color='grey', linewidth=2.0)

    # We change the fontsize of minor ticks label
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.tick_params(axis='both', which='minor', labelsize=8)

    # Plot the usage of ECMs as vectors
    for i in range(len(Xusage)):
        color = cmap_curr[i, :]
        ax.quiver(Xusage[i], Yusage[i], Uusage[i], Vusage[i], pivot='tail', angles='xy', scale_units='xy',
                  linestyle='--', color=color, scale=1, width=quiverwidth)

    # if objective_dict['obj_cons'] == 'obj':
    #     ax.set_title("Fraction needed \n per %.2f: " % scaling_factor + obj_string, size=8)
    #     # plt.title("Fraction needed per %.2f:\n" % scaling_factor + obj_string)
    # else:
    #     ax.set_title("Fraction needed for \n demanded: " + obj_string, size=8)
    #     # plt.title("Fraction needed for demanded:\n" + obj_string)

    # custum legend: plot active ECMs for which also a vector is plotted.
    legend_elements = []
    for i in curr_act_inds:
        # Only plot if there is a vector with non-zero length for the current set of constraints
        if Uusage[i]!=0 or Vusage[i]!=0:
            # lines
            #legend_elements.append(Line2D([0], [0], color=cmap_curr[i, :], #'w', markerfacecolor=cmap_curr[i, :], marker='o',
            #                          label=actECMs.index[i])) #, markersize=scatter_normal/3*2, alpha=alpha_active))
            # scatter dots
            legend_elements.append(
                Line2D([0], [0], color='w', markerfacecolor=cmap_curr[i, :], marker='o',
                       label=actECMs.index[i], markersize=scatter_normal/3*2, alpha=alpha_active))

    ax.legend(handles=legend_elements, bbox_to_anchor=(1., 0.5), loc='center left', prop={'size': 8},
              title='ECM ID', frameon=False,)
    #ax.legend(loc='upper right')
    # set legend of axes invisible
    # ax.get_legend().set_visible(False)
    # get legend handles and labels to create shared legend
    handles, labels = ax.get_legend_handles_labels()

    # Create common x and y axis labels
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)

    plt.xlabel('Needed fraction of constraint: ' + x_label)
    if not 'virtual' in cons_IDs:
        plt.ylabel(ylabel='Needed fraction of constraint: ' + y_label)

    # Create shared legend based on latest ax handles and labels
    fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.), loc='lower center',
               borderaxespad=0., ncol=2,
               frameon=False, prop={'size': 8})
    plt.subplots_adjust(bottom=0.15, hspace=0.6)

    if result_dir:
        fig.savefig(
            os.path.join(result_dir,
                         "cost_plot_ECMs" + model_dict['model_name'] + "_" + cons_IDs[0] + "_" + cons_IDs[1] + ".png"),
            bbox_inches="tight")
    else:
        fig.savefig(
            os.path.join(os.getcwd(),
                         "cost_plot_ECMs" + model_dict['model_name'] + "_" + cons_IDs[0] + "_" + cons_IDs[1] + ".png"),
            bbox_inches="tight")

