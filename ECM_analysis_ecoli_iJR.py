import pickle
import itertools

from helpers_EFM_FBA import *

"""CONSTANTS"""
# TODO: Try a bigger model
# model_name = "e_coli_core"
# model_name = "iAB_RBC_283"
# model_name = "iNF517" #Flahaut lactis model Bigg version
model_name = "iJR904"
# model_name = "iAF1260b"
#model_name = "iAF1260b"  # E.coli
#model_name = "iAF1260b_new" # E.coli
#model_name = "lplawcfs1" #plantarum model Bas # species ending with _b are boundary species. extracellular compartment
#model_name = "Lactobacillus_plantarum_WCFS1_Official_23_May_2019_18_45_01"
# denoted with 'Extra_organism'
#model_name = "MG1363_20190628"  # lactis MG1363 model as updated for pcLactis, ECMs can be calculated for active, activities not retrieved
# network, but something goes wrong when calculating the activities of the ECMs in the FBA solution. Supremum norm is non-zero.
#model_name = "iIT341" # Helicobacter pylori 26695 works.
#model_name = "iYO844" # Bacillus subtilis subsp. subtilis str. 168; works

# Define directories for storing: 1) newly made models, 2) results such as figures, and 3) for finding models
sbml_dir = os.path.join(os.getcwd(), "data", model_name + "_2", 'models', 'sbml')
result_dir = os.path.join(os.getcwd(), "data", model_name + "_2", "EFM_yield_analysis")
model_dir = os.path.join("models")
LOAD_IF_AVAILABLE = False  # If this is true, newly created models will be loaded if the program runs a second time
N_KOS = 0  # number of knockout models that will be created
# adjust DROP_MODELS based on this.
# The complete model is sometimes to large to calculate ECMs. This is therefore dropped manually
DROP_MODEL_TAGS = ['full', 'active', 'fva', 'active_hidden', 'ko']  # ['full','active','fva','hidden','active_hidden' ,'ko']
USE_EXTERNAL_COMPARTMENT = None
ONLY_RAYS = False  # True or False
SAVE_result = False  # saves list_model_dict after ECM enumeration and calculation of ECM activities in FBA solution
VERBOSE = True

"""START ANALYSIS"""
"""Create folders for storing the results"""
if not os.path.exists(sbml_dir):
    os.makedirs(sbml_dir)
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

"""Load main model and clone it to have a model that we can adapt"""
try:
    cmod = cbm.readSBML3FBC(os.path.join(model_dir, model_name + ".xml"))
except:  # use version 2
    cmod = cbm.readSBML2FBA(os.path.join(model_dir, model_name + ".xml"))

cbm.doFBAMinSum(cmod)
intermediate_cmod = cmod.clone()

"""Adapt model to your liking, and update the FBA"""
oxygen_reaction = intermediate_cmod.getReaction('R_EX_o2_e')
oxygen_reaction.setLowerBound(-10)  # -15, -10, -5-> give two active constraints; -1000 -> give 1 active constraint
atp_reaction = intermediate_cmod.getReaction('R_ATPM')
atp_reaction.setLowerBound(0.)
intermediate_cmod.setReactionLowerBound("R_EX_co2_e", 0.)
intermediate_cmod.setReactionLowerBound("R_EX_fe2_e", 0.)
intermediate_cmod.setReactionLowerBound("R_EX_h2o_e", 0.)
intermediate_cmod.setReactionLowerBound("R_EX_h_e", 0.)
intermediate_cmod.setReactionLowerBound("R_EX_k_e", 0.)
intermediate_cmod.setReactionLowerBound("R_EX_na1_e", 0.)

# Find FBA objective value
original_FBA_val = cbm.doFBA(intermediate_cmod)

"""Filter infeasible reactions in intermediate_cmod (reactions with zero lower and upper bound)."""
species = list(intermediate_cmod.species)
species_index = {item.id: index for index, item in enumerate(species)}
reactions = intermediate_cmod.reactions

# Check if all reactions can run
not_feasible_list = []
for reaction in reactions:
    reaction_feasible = True
    lowerBound, upperBound, _ = intermediate_cmod.getFluxBoundsByReactionID(reaction.id)
    if lowerBound is not None and upperBound is not None:
        if lowerBound.value == 0 and upperBound.value == 0:
            reaction_feasible = False
        elif lowerBound.value > upperBound.value:
            reaction_feasible = False

    elif lowerBound is None and upperBound is None: # check the equality notion at the end of getReactionBounds
        # When saved as sbml3FBC2 these bounds will be saved as 0.0 bounds.
        reaction_feasible = False

    if not reaction_feasible:
        not_feasible_list.append(reaction.id)
        print('Reaction {} has zero lower bound and upper bound or equal to None, '
              'and is therefore infeasible. Please adjust bounds or delete reaction.'.format(
            reaction.id))
        # raise Exception(
        #    'Reaction {} has lower bound {} and upper bound {}, '
        #    'and is therefore infeasible. Please adjust bounds or delete reaction.'.format(
        #        reaction.id, lowerBound.value, upperBound.value))

for rid in not_feasible_list:
    intermediate_cmod.deleteReactionAndBounds(rid)

"""Check reversibility of reactions. If lb<0 than reaction should be reversible for ecmtool. 
This is according to SBML language rules."""
for rid in intermediate_cmod.getReactionIds():
    if intermediate_cmod.getReactionLowerBound(rid) < 0 and not intermediate_cmod.getReaction(rid).reversible:
        print('Reversibility of ' + rid + ' is set to True because lower bound was negative.')
        intermediate_cmod.getReaction(rid).reversible = True

cbm.doFBA(intermediate_cmod)
cbm.doFBAMinSum(intermediate_cmod)

"""Create list with a dictionary in which all info will be stored for each model"""
list_model_dicts = []

# Create dictionary for original model
list_model_dicts.append(
    {'model': intermediate_cmod, 'model_name': 'original_network', 'calc_efms': False, 'get_activities': True,
     'get_relevant_efms': False, 'drop_tag': 'full'})

# Deletes all non active reactions (determined by seeing if the flux is lower than a tolerance)
sub_model_name = 'active_subnetwork'
active_model_path = os.path.join(sbml_dir, sub_model_name + ".xml")
if os.path.exists(active_model_path) and LOAD_IF_AVAILABLE:  # Checks if this model was already made
    try:
        intermediate_cmod_active = cbm.readSBML3FBC(active_model_path)
        # intermediate_cmod_active.setNotes(str(intermediate_cmod_active.getNotes().encode(encoding='ascii', errors='ignore')))
    except:  # use version 2
        intermediate_cmod_active = cbm.readSBML2FBA(active_model_path)
        # intermediate_cmod_active.setNotes(str(intermediate_cmod_active.getNotes().encode(encoding='ascii', errors='ignore')))
else:
    intermediate_cmod_active = delete_non_active_network(intermediate_cmod, which_zeros='flux_zero', zero_tol=1e-15,
                                                         opt_tol=1e-8)  # 15, 8
    # intermediate_cmod_active.setNotes(str(intermediate_cmod_active.getNotes().encode(encoding='ascii', errors='ignore')))

# Create dictionary for active subnetwork
list_model_dicts.append(
    {'model': intermediate_cmod_active, 'model_name': sub_model_name, 'calc_efms': False, 'get_activities': True,
     'get_relevant_efms': False, 'drop_tag': 'active'})

# Delete only reactions that are never active (based on FVA)
if 'fva' not in DROP_MODEL_TAGS:
    sub_model_name = 'active_subnetwork_FVA'
    model_path = os.path.join(sbml_dir, sub_model_name + ".xml")
    if os.path.exists(model_path) and LOAD_IF_AVAILABLE:  # Checks if this model was already made
        try:
            intermediate_cmod_FVA_active = cbm.readSBML3FBC(model_path)
            intermediate_cmod_FVA_active.setNotes(
                intermediate_cmod_FVA_active.getNotes().encode(encoding='ascii', errors='ignore'))
        except:  # use version 2
            intermediate_cmod_FVA_active = cbm.readSBML2FBA(model_path)
    else:
        intermediate_cmod_FVA_active = delete_non_active_network(intermediate_cmod, which_zeros='FVA_zero')

    # Create dictionary for active subnetwork based on FVA
    list_model_dicts.append(
        {'model': intermediate_cmod_FVA_active, 'model_name': sub_model_name, 'calc_efms': False,
         'get_activities': True, 'get_relevant_efms': False,  'drop_tag': 'fva'})

# Optional: Create list of models with reactions that are active when certain knockouts are created
# This is a way to create more models to try-out stuff. You can use this option by setting n_KOs > 0
list_KO_models = create_KO_models(intermediate_cmod, n_KOs=N_KOS)
for model_ind, model in enumerate(list_KO_models):
    list_model_dicts.append({'model': model, 'model_name': 'active_subnetwork_KO_%d' % model_ind, 'calc_efms': True,
                             'get_activities': True, 'get_relevant_efms': True,  'drop_tag': 'ko'}) #was False

"""Rebuild stoichiometric matrix for reduced models"""
for model_dict in list_model_dicts:
    model_dict['model'].buildStoichMatrix()
    cbm.doFBAMinSum(model_dict['model'])

"""Store the created models in the result- and in the xml-directory"""
for model_dict in list_model_dicts:
    model_id = model_dict['model_name']
    model = model_dict['model']
    model_path = os.path.join(sbml_dir, model_id + ".xml")
    model_dict['model_path'] = model_path
    if not os.path.exists(model_path) or not LOAD_IF_AVAILABLE:
        cbm.writeSBML3FBCV2(model, model_path, add_cbmpy_annot=False)
        cbm.writeSBML3FBCV2(model, os.path.join(result_dir, model_id + ".xml"), add_cbmpy_annot=False)

full_model_path = list_model_dicts[0]['model_path']

"""Find the relevant, i.e., growth-limiting, constraints for the original model"""
# Do non-zero reduced cost analysis
cbm.analyzeModel(intermediate_cmod)
nzrc_dictionaries, n_objectives = get_nzrc(intermediate_cmod)
cbm.doFBAMinSum(intermediate_cmod)

"""Check which metabolite is coupled to the constrained reactions. If no metab is coupled, tag the reaction
with a virtual metabolite"""
nzrc_dictionaries, reactions_to_tag = findConstraintMetabolites(nzrc_dictionaries, intermediate_cmod)

"""Determine active objective function. Split the non-zero reduced costs in objectives and constraints"""
infos_obj, infos_cons = get_info_objectives_constraints(nzrc_dictionaries, intermediate_cmod)

"""Get relevant metabolites for the cost-calculations, find indices of external metabs that can be ignored"""
# We only have to focus on conversions of metabolites that are somehow active in a constraint
relevant_metabs = [info_dict['ext_metab'] for info_dict in infos_obj + infos_cons]
# relevant_metabs = ['M_glc__D_e', 'M_o2_e', 'objective']
# The following finds all indices of external metabolites that can be ignored in the conversion analysis
hide_indices = find_hide_indices(full_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs,
                                 use_external_compartment=USE_EXTERNAL_COMPARTMENT)
# By default no metabolites are hidden
for model_dict in list_model_dicts:
    model_dict['hide_metabolites'] = []

# try out with only glucose and objective
# relevant_metabs = ['objective', 'M_glc__D_e']
# hide_indices = find_hide_indices(full_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs,
#                                     use_external_compartment = USE_EXTERNAL_COMPARTMENT)
"""Create another dictionary for a model in which some metabolites are hidden"""
list_model_dicts.append(
    {'model': intermediate_cmod, 'model_name': 'original_with_hidden_metabolites', 'calc_efms': False,
     'get_activities': True, 'hide_metabolites': hide_indices, 'get_relevant_efms': True,
     'model_path': os.path.join(sbml_dir, "original_network.xml"), 'drop_tag': 'hidden'})

"""Create another dictionary for the active model with some metabolites hidden"""
active_hide_indices = find_hide_indices(active_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs,
                                        use_external_compartment=USE_EXTERNAL_COMPARTMENT)
# TODO: adjust below !!! modelpath etc...............................................................
list_model_dicts.append(
    {'model': intermediate_cmod_active, 'model_name': 'active_network_with_hidden_metabolites', 'calc_efms': False, #False
     'get_activities': True, 'hide_metabolites': active_hide_indices, 'get_relevant_efms': True,
     'model_path': active_model_path,  'drop_tag': 'active_hidden'})


"""Calculate ECMs and/or EFMs"""
# For genome-scale models we cannot yet calculate ECMs, then we should calculate them for only the active subnetwork,
# or only for the model with external metabolites that are ignored.

# We remove some models defined at top of the file
list_model_dicts_remember = list_model_dicts.copy()
list_model_dicts = [model_dict for model_dict in list_model_dicts if model_dict['drop_tag'] not in DROP_MODEL_TAGS]

# Load ECMs file E.coli groot iJR
#ecms_df_pre = pd.read_csv(os.path.join(result_dir, 'ECMs_hiddennetwork.csv')) # rounded to two decimals
ecms_df_pre = pd.read_csv(os.path.join(result_dir, 'iJR_hideallexceptglco2biomass.csv')) # not rounded
ecms_df_pre = ecms_df_pre.transpose()

# and update model_dict
ecms_from_efms = None
efms_df = None
#full_network_ecm = None
for model_dict in list_model_dicts:
    full_network_ecm = get_full_network(file_path=model_dict['model_path'],
                                        reactions_to_tag=reactions_to_tag,
                                        print_results=False,
                                        hide_metabs=model_dict['hide_metabolites'],
                                        use_external_compartment=USE_EXTERNAL_COMPARTMENT,
                                        only_rays=ONLY_RAYS
                                        )
    # Find all metabolite_ids in the order used by ECMtool
    mets = [met.id for met in full_network_ecm.metabolites]

    ecms_df_pre_norm = pd.DataFrame(np.zeros((len(mets), len(ecms_df_pre.columns))), index=mets, columns=ecms_df_pre.columns)
    for met in ecms_df_pre.index:
        ecms_df_pre_norm.loc[met] = ecms_df_pre.loc[met]

    # Normalize the ECMs. We first normalize all ECMs that produce the objective to produce one objective, after that
    # we normalize to different metabolite productions, such that the ECMs are maximally comparable
    ecms_matrix = normalize_ECMS_objective_first(ecms_df_pre_norm.to_numpy(), full_network_ecm, infos_obj)

    # Create dataframe with the ecms as columns and the metabolites as index
    ecms_df = pd.DataFrame(ecms_matrix, index=mets)

    model_dict.update(
        {'ecms_df': ecms_df, 'network': full_network_ecm, 'efms_df': efms_df, 'ecms_from_efms': ecms_from_efms})

"""Select only the relevant part of the ECM information. So, only the consumption/production flux of metabolites that
    are associated to a flux bound."""
for model_dict in list_model_dicts:
    #filtered_infos_cons = [rid for rid in infos_cons if rid['rid'] in model_dict['model'].getReactionIds()]
    filtered_infos_cons = infos_cons
    filtered_infos_obj = infos_obj
    #filtered_infos_obj = [rid for rid in infos_obj if rid['rid'] in model_dict['model'].getReactionIds()]
    model_dict['infos_cons'] = filtered_infos_cons
    model_dict['infos_obj'] = filtered_infos_obj
    #print('infos_cons and infos_obj were filtered')
    #table_cons_df, table_obj_df = get_relevant_parts_ECMs(model_dict['ecms_df'], infos_cons, infos_obj)
    table_cons_df, table_obj_df = get_relevant_parts_ECMs(model_dict['ecms_df'], model_dict['infos_cons'], model_dict['infos_obj'])


    # Then calculate the activities of the ECMs in the FBA solution and add it to the dataframes.
    if model_dict['get_activities']:
        print('Calculating activities for ' + model_dict['model_name'])
        activities = get_activity_ECMs(model_dict['ecms_df'], model_dict['model'],
                                       model_dict['network'], hide_indices=model_dict['hide_metabolites'])
        # Store the activities of the ECMs in the tables where we had the relevant info
        table_cons_df['active'] = activities
        table_obj_df['active'] = activities

    model_dict['table_cons_df'] = table_cons_df
    model_dict['table_obj_df'] = table_obj_df

"""Plot these costs"""
# We need the objective value to scale the axes of our plots
objective_val = intermediate_cmod.getObjFuncValue()
# We make 2-D cost plots only (3-D is too confusing). Since each constrained metabolite brings costs, each metabolite
# needs to be plotted on some axes. We thus make ceil(n_cons/2) figures
constrained_ext_metabs = [cons_dict['ext_metab'] for cons_dict in list_model_dicts[0]['infos_cons']]
if len(constrained_ext_metabs) > 1:
    if model_name == "e_coli_core":
        # all possible combinations of constraints
        cons_ID_combi_list = list(itertools.combinations(constrained_ext_metabs, 2))
    else:
        cons_ID_combi_list = get_cons_ID_combinations(constrained_ext_metabs, list_model_dicts[0])
else:
    cons_ID_combi_list = None

if 'MG1363' in model_name:
    cons_ID_combi_list += [['M_asn__L_e', 'M_asp__L_e'], ['M_glc__D_e', 'M_ile__L_e']] # ['M_asp__L_e', 'M_asn__L_e']

if cons_ID_combi_list:
    for cons_ID_combi in cons_ID_combi_list:
        # Plot the results of the different models (with different (sub)networks) in one plot
        plot_different_approaches_one_figure(list_model_dicts, infos_obj,  # list_model_dicts_remember[0]['infos_obj'],
                                             infos_cons,  # list_model_dicts_remember[0]['infos_cons'],
                                             obj_val=objective_val,
                                             result_dir=result_dir, cons_IDs=cons_ID_combi)
else:
    #for ind in range(0, len(constrained_ext_metabs), 2):
    # print(ind)
    # Select constraints that will be plotted in this figure
    #cons_for_fig = constrained_ext_metabs[ind:max(ind + 2, len(constrained_ext_metabs))]
    cons_for_fig = constrained_ext_metabs
    # Plot the results of the different models (with different (sub)networks) in one plot
    plot_different_approaches_one_figure(list_model_dicts, infos_obj,  # list_model_dicts_remember[0]['infos_obj'],
                                             infos_cons,  # list_model_dicts_remember[0]['infos_cons'],
                                             obj_val=objective_val,
                                             result_dir=result_dir, cons_IDs=cons_for_fig)

# Now make for each approach for which we calculated the activities of the ECMs a different plot in which the usage of
# the ECMs is shown with vectors.
if cons_ID_combi_list:
    for model_dict in list_model_dicts:
        if model_dict['get_activities']:
            print('Plotting the cost vectors including usage for model %s' % model_dict['model_name'])
            for cons_ID_combi in cons_ID_combi_list:
                print(cons_ID_combi)
                plot_costs(model_dict, infos_obj, infos_cons,  # model_dict['infos_obj'], model_dict['infos_cons'],
                           cons_IDs=cons_ID_combi, obj_val=objective_val,
                           show_active=True, result_dir=result_dir)
                plot_costs_ECMs(model_dict, infos_cons, cons_IDs=cons_ID_combi, result_dir=result_dir)
else:
    for model_dict in list_model_dicts:
        if model_dict['get_activities']:
            print('Plotting the cost vectors including usage for model %s' % model_dict['model_name'])
            plot_costs(model_dict, infos_obj, infos_cons,  # model_dict['infos_obj'], model_dict['infos_cons'],
                       cons_IDs=constrained_ext_metabs, obj_val=objective_val,
                       show_active=True, result_dir=result_dir)
            plot_costs_ECMs(model_dict, infos_cons, cons_IDs=constrained_ext_metabs, result_dir=result_dir)

"""Find corresponding EFM(s) for active ECM(s)"""
# TODO: Make a function that finds an EFM corresponding to an ECM
# TODO: Apply this function to the resulting ECMs for the hidden metabolites, to get a minimal network needed for these constraints
# TODO: Apply this function (if feasible) for all ECMs, to get a minimal FBA-network

# Find some EFM that corresponds to each of the active ECMs in the hidden-network
for model_dict in list_model_dicts:
    if model_dict['get_relevant_efms'] and model_dict['get_activities']:
        relevant_efms_df, full_relevant_ecms_df = find_associated_efms(model_dict['model'], model_dict['table_cons_df'],
                                                                       model_dict['ecms_df'],
                                                                       infos_obj + infos_cons, model_dict['model_path'],
                                                                       use_external_compartment=USE_EXTERNAL_COMPARTMENT,
                                                                       only_relevant_ECMs=False)
        relevant_efms_df.to_csv(
            os.path.join(
                result_dir, "efms_corresponding_to_hide_ecms" + model_dict['model_name'] + ".csv")) #, index=False)
        full_relevant_ecms_df.to_csv(
            os.path.join(
                result_dir, "full_ecms_corresponding_to_hide_ecms" + model_dict['model_name'] + ".csv")) #, index=False)

        if VERBOSE:
            print_ecms_direct(np.transpose(full_relevant_ecms_df.values), full_relevant_ecms_df.columns)


# Todo: find differences/similarities between ECMs and EFMs
a = full_relevant_ecms_df.iloc[0]-full_relevant_ecms_df.iloc[1]
a.to_numpy().nonzero()
full_relevant_ecms_df.iloc[:,a.to_numpy().nonzero()[0]]

b = relevant_efms_df.iloc[0] - relevant_efms_df.iloc[1]
b.to_numpy().nonzero()
relevant_efms_df.iloc[:,b.to_numpy().nonzero()[0]]

print(relevant_efms_df.iloc[:,b.to_numpy()>1.].transpose())

for rid in relevant_efms_df.iloc[:,b.to_numpy()>1.].transpose().index:
    reaction = intermediate_cmod.getReaction(rid)
    print(reaction.getName())
    print(relevant_efms_df[rid])
    print(reaction.getEquation())

#
result_nzrc = pd.DataFrame()
for nzrc in nzrc_dictionaries:
    print(nzrc["rid"],nzrc["nzrc"], nzrc["flux_val"])
    result_nzrc[nzrc["rid"]] = [nzrc["nzrc"], nzrc["flux_val"]]
result_nzrc.index = ['nzrc', 'flux_val']
#result_nzrc.transpose()

# TODO: plot the EFMs on a map of the metabolic network - or overview/lumped metabolic network
# Idea: select reactions for the map based on the active reactions in the EFMs.
# Lump reactions and model if needed for visualization. <- still to be developed
# Create groups of reactions based on EFMs, if in linear pathway add a lumped reaction to model and to group and remove
# reactions that are replaced with the lumped variant.
# Color based on activity (in case of 1 EFM)
# or heatmap of EFMs with activity?



"""Problem solving"""
for model_dict in list_model_dicts:
    print('Model =', model_dict['model_name'])
    print('Activity per ECM for objective(s)')
    print(model_dict['table_obj_df'][model_dict['table_obj_df']['active']!=0.])
    print('Activity per ECM for the active constraint(s)')
    print(model_dict['table_cons_df'][model_dict['table_cons_df']['active']!=0.])
    #print('Glucose flux:', sum(model_dict['table_cons_df']['active']*model_dict['table_cons_df']['M_glc__D_e']))
    #print('Oxygen flux:', sum(model_dict['table_cons_df']['active']*model_dict['table_cons_df']['M_o2_e']))
    print()
    print('Active ECMs and used metabolites')
    print(get_active_ecms(model_dict))
    print()