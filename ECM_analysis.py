import pickle

from helpers_EFM_FBA import *

"""CONSTANTS"""
# TODO: Try a bigger model
model_name = "e_coli_core"
# With the following models, some shit goes wrong always
# model_name = "iAB_RBC_283"
# model_name = "iNF517" #Flahaut lactis model Bigg version
# model_name = "iJR904"
# model_name = "iAF1260b"
#model_name = "iAF1260b"  # E.coli
#model_name = "iAF1260b_new" # E.coli
#model_name = "lplawcfs1" #plantarum model Bas # species ending with _b are boundary species. extracellular compartment
#model_name = "Lactobacillus_plantarum_WCFS1_Official_23_May_2019_18_45_01"
# denoted with 'Extra_organism'
model_name = "MG1363_20190628"  # lactis MG1363 model as updated for pcLactis, ECMs can be calculated for active, activities not retrieved
# network, but something goes wrong when calculating the activities of the ECMs in the FBA solution. Supremum norm is non-zero.
#model_name = "iIT341" # Helicobacter pylori 26695 works.
#model_name = "iYO844" # Bacillus subtilis subsp. subtilis str. 168; works

# Define directories for storing: 1) newly made models, 2) results such as figures, and 3) for finding models
sbml_dir = os.path.join(os.getcwd(), "data", model_name + "_2", 'models', 'sbml')
result_dir = os.path.join(os.getcwd(), "data", model_name + "_2", "EFM_yield_analysis")
model_dir = os.path.join("models")
LOAD_IF_AVAILABLE = False  # If this is true, newly created models will be loaded if the program runs a second time
N_KOS = 0 # number of knockout models that will be created
# adjust DROP_MODELS based on this.
#DROP_MODELS = [0,12]  # The complete model is sometimes to large to calculate ECMs. This is therefore dropped manually
DROP_MODELS = [0, 2]  # The complete model is sometimes to large to calculate ECMs. This is therefore dropped manually
EXTERNAL_COMPARTMENT = 'e'

"""START ANALYSIS"""
"""Create folders for storing the results"""
if not os.path.exists(sbml_dir):
    os.makedirs(sbml_dir)
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

"""Load main model and clone it to have a model that we can adapt"""
try:
    cmod = cbm.readSBML3FBC(os.path.join(model_dir, model_name + ".xml"))
    # cmod.setNotes(cmod.getNotes().encode(encoding='ascii', errors='ignore'))
except:  # use version 2
    cmod = cbm.readSBML2FBA(os.path.join(model_dir, model_name + ".xml"))

cbm.doFBAMinSum(cmod)
intermediate_cmod = cmod.clone()
# intermediate_cmod.setNotes(str(intermediate_cmod.getNotes()))


"""Adapt model to your liking, and update the FBA"""
if model_name == "e_coli_core":
    # E. Coli Core specific adaptations
    # Set own constraints (no constraint on ATP, constraint on O2 though)
    atp_reaction = intermediate_cmod.getReaction('R_ATPM')
    # atp_reaction.setLowerBound(0.) ## eunice_edit 0.
    oxygen_reaction = intermediate_cmod.getReaction('R_EX_o2_e')
    oxygen_reaction.setLowerBound(-15.)
    acetate_reaction = intermediate_cmod.getReaction('R_EX_ac_e')
    #acetate_reaction.setLowerBound(-5.) # This works with ub=1000; reaction set reversible.
    #acetate_reaction.setUpperBound(-1.)
    #acetate_reaction.setUpperBound(0.) # This works with new version network.py with lb=-5.; reaction set reversible.
    #acetate_reaction.setLowerBound(2.) # This works with ub=5; reaction set not reversible
    #acetate_reaction.setUpperBound(5.) # This works with lb=2; reaction set not reversible
    acetate_reaction.reversible=True
    glc_reaction = intermediate_cmod.getReaction('R_EX_glc__D_e')
    glc_reaction.setUpperBound(-8.)
    glc_reaction.setLowerBound(-10.) #is original
    #glc_reaction.reversible = False # was True

    DROP_MODELS = [0, 2] #[0, 2]  # [0]  # The complete model is sometimes to large to calculate ECMs. This is therefore dropped manually

elif model_name == 'iAF1260b_new':
    # E. coli
    oxygen_reaction = intermediate_cmod.getReaction('R_EX_o2_e')
    oxygen_reaction.setLowerBound(-10.)  # -> leads to two active constraints
    atp_reaction = intermediate_cmod.getReaction('R_ATPM')
    #atp_reaction.setLowerBound(.)

elif model_name == 'MG1363_20190628':
    GOEL_BOUNDS = 'basic'  # either 'basic', 'FBA', 'FBA_adjusted'
    # L. lactis MG1363
    # Read excel file with measured exchange bounds
    filepath = os.path.join(result_dir, "..", "goel_bounds.xlsx")
    goel_bounds = pd.read_excel(filepath, index_col=1, sheet_name='For FBA')
    goel_bounds_adjusted = pd.read_excel(filepath, index_col=1, sheet_name='For FBA adjusted')

    # Set exchange reaction bounds to basis rids
    filepath = os.path.join(result_dir, "..", "exchange_rids_basic.txt")
    exchange_rids_basic = pd.read_csv(filepath, sep='\t', index_col=0)
    exchange_rids_basic = exchange_rids_basic.astype({'lb': float, 'ub': float})
    if GOEL_BOUNDS == 'basic':
        for rid in exchange_rids_basic.index:
            # print(cmod.getReactionBounds(rid))
            intermediate_cmod.setReactionBounds(rid, exchange_rids_basic["lb"][rid], exchange_rids_basic["ub"][rid])
            # print(intermediate_cmod.getReactionBounds(rid))

        # if using basis bounds, also set glucose lower bound open.
        intermediate_cmod.setReactionBounds("R_EX_glc__D_e", -4., 0.)
    elif GOEL_BOUNDS == 'FBA':
        # Set exchange bounds to dilution rate of choice
        # This stil doesn't work. ecmtool errors, and takes a lot of time... At least several hours.
        for rid in goel_bounds.index:
            if rid not in []:#["R_EX_pi_e", "R_EX_nac_e", "R_EX_o2_e", "R_EX_orn_e"]:#, "R_EX_pyr_e", "R_EX_succ_e", "R_EX_nh4_e", "R_EX_lac__L_e"]:
                # 0.6
                print(intermediate_cmod.getReactionBounds(rid))
                intermediate_cmod.setReactionBounds(rid, goel_bounds['LB0.6A'][rid], goel_bounds['UB0.6A'][rid])
                print(intermediate_cmod.getReactionBounds(rid))
    elif GOEL_BOUNDS == 'FBA_adjusted':
        for rid in goel_bounds_adjusted.index:
            if rid not in []:#["R_EX_pi_e", "R_EX_nac_e", "R_EX_o2_e", "R_EX_orn_e"]:#, "R_EX_pyr_e", "R_EX_succ_e", "R_EX_nh4_e", "R_EX_lac__L_e"]:
                # 0.6
                print(intermediate_cmod.getReactionBounds(rid))
                intermediate_cmod.setReactionBounds(rid, goel_bounds_adjusted['LB0.6A'][rid], goel_bounds_adjusted['UB0.6A'][rid])
                print(intermediate_cmod.getReactionBounds(rid))
        intermediate_cmod.setReactionUpperBound('R_EX_pnto__R_e', 0.) # make input
        intermediate_cmod.setReactionUpperBound('R_EX_thm_e', 0.) #make input
        for rid in ['R_EX_h2o_e', 'R_EX_h2s_e', 'R_EX_h_e']:
            # make output
            intermediate_cmod.setReactionLowerBound(rid, 0.)
            #

    DROP_MODELS = [0, 2, 3]

elif model_name == "lplawcfs1": #plantarum model Bas
    # something goes wrong here. biomass reaction has an equal constraint. you cannot set this with python 3 cbmpy.
    # This causes errors later on. upper and lower bound are None, equal is 0.0, then reaction is 'infeasible'.
    # you cannot change the equal bound or remove it to set a lower or upper bound.
    intermediate_cmod.createObjectiveFunction('R_biomass_LPL60')
    EXTERNAL_COMPARTMENT = 'Extern_organism'

elif model_name == "Lactobacillus_plantarum_WCFS1_Official_23_May_2019_18_45_01":
    for rid in intermediate_cmod.getReactionIds(): #"R_EX_"):
        print(intermediate_cmod.getReactionBounds(rid))
        intermediate_cmod.setReactionBounds(rid, intermediate_cmod.getReactionLowerBound(rid) * 100.,
                                                 intermediate_cmod.getReactionUpperBound(rid) * 100.)

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
        print(rid)
        intermediate_cmod.getReaction(rid).reversible = True

cbm.doFBA(intermediate_cmod)
cbm.doFBAMinSum(intermediate_cmod)

"""Create list with a dictionary in which all info will be stored for each model"""
list_model_dicts = []

# Create dictionary for original model
list_model_dicts.append(
    {'model': intermediate_cmod, 'model_name': 'original_network', 'calc_efms': False, 'get_activities': True,
     'get_relevant_efms': False})

# Deletes all non active reactions (determined by seeing if the flux is lower than a tolerance)
model_name = 'active_subnetwork'
model_path = os.path.join(sbml_dir, model_name + ".xml")
if os.path.exists(model_path) and LOAD_IF_AVAILABLE:  # Checks if this model was already made
    try:
        intermediate_cmod_active = cbm.readSBML3FBC(model_path)
        # intermediate_cmod_active.setNotes(str(intermediate_cmod_active.getNotes().encode(encoding='ascii', errors='ignore')))
    except:  # use version 2
        intermediate_cmod_active = cbm.readSBML2FBA(model_path)
        # intermediate_cmod_active.setNotes(str(intermediate_cmod_active.getNotes().encode(encoding='ascii', errors='ignore')))
else:
    intermediate_cmod_active = delete_non_active_network(intermediate_cmod, which_zeros='flux_zero', zero_tol=1e-15,
                                                         opt_tol=1e-8)  # 15, 8
    # intermediate_cmod_active.setNotes(str(intermediate_cmod_active.getNotes().encode(encoding='ascii', errors='ignore')))

# Create dictionary for active subnetwork
list_model_dicts.append(
    {'model': intermediate_cmod_active, 'model_name': model_name, 'calc_efms': True, 'get_activities': True,
     'get_relevant_efms': True})

# Delete only reactions that are never active (based on FVA)
model_name = 'active_subnetwork_FVA'
model_path = os.path.join(sbml_dir, model_name + ".xml")
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
    {'model': intermediate_cmod_FVA_active, 'model_name': model_name, 'calc_efms': False,
     'get_activities': True, 'get_relevant_efms': False})

# Optional: Create list of models with reactions that are active when certain knockouts are created
# This is a way to create more models to try-out stuff. You can use this option by setting n_KOs > 0
list_KO_models = create_KO_models(intermediate_cmod, n_KOs=N_KOS)
for model_ind, model in enumerate(list_KO_models):
    list_model_dicts.append({'model': model, 'model_name': 'active_subnetwork_KO_%d' % model_ind, 'calc_efms': False,
                             'get_activities': True, 'get_relevant_efms': False}) #was False

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
# The following finds all indices of external metabolites that can be ignored in the conversion analysis
hide_indices = find_hide_indices(full_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs)
# By default no metabolites are hidden
for model_dict in list_model_dicts:
    model_dict['hide_metabolites'] = []

# try out with only glucose and objective
#relevant_metabs = ['objective', 'M_glc__D_e']
#hide_indices = find_hide_indices(full_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs)
"""Create another dictionary for a model in which some metabolites are hidden"""
list_model_dicts.append(
    {'model': intermediate_cmod, 'model_name': 'original_with_hidden_metabolites', 'calc_efms': False,
     'get_activities': True, 'hide_metabolites': hide_indices, 'get_relevant_efms': False,
     'model_path': os.path.join(sbml_dir, "original_network.xml")})

"""Create another dictionary for the active model with some metabolites hidden"""
if model_name == '':
    active_model_path = list_model_dicts[0]['model_path']
    hide_indices = find_hide_indices(active_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs)
    # TODO: adjust below !!! modelpath etc...............................................................
    list_model_dicts.append(
        {'model': intermediate_cmod_active, 'model_name': 'active_network_with_hidden_metabolites', 'calc_efms': False,
         'get_activities': True, 'hide_metabolites': hide_indices, 'get_relevant_efms': False,
         'model_path': os.path.join(sbml_dir, "original_network.xml")})

# Debug?! get network and print external metabolites.
# file_path = os.path.join(sbml_dir, 'active_subnetwork' + ".xml")
# network = get_network_class(file_path)
# external_met_IDs = [network.metabolites[i].id for i in network.external_metabolite_indices()]
# for i in network.external_metabolite_indices():
#     print(network.metabolites[i].id, network.metabolites[i].direction)
#     print(intermediate_cmod_active.getSpecies(network.metabolites[i].id).isReagentOf())
#     for rid in intermediate_cmod_active.getSpecies(network.metabolites[i].id).isReagentOf():
#         if "_EX_" in rid:
#             print(intermediate_cmod_active.getReaction(rid).value)
# print(external_met_IDs)
# external_model_sids = [species.id for species in intermediate_cmod_active.species if species.compartment == 'e']
# print(external_model_sids)
#
# print(set(external_met_IDs) - set(external_model_sids))
# print(set(external_model_sids) - set(external_met_IDs))

#print(set(external_met_IDs) - set(basic_external_met_IDs))
#print(set(basic_external_met_IDs) - set(external_met_IDs))
#basic_external_met_IDs = external_met_IDs
#basic_external_model_sids = external_model_sids

#file_path_original = os.path.join(sbml_dir, 'original_network' + ".xml")
#network_original = get_network_class(file_path_original)
#external_met_IDs_original = [network_original.metabolites[i].id for i in network_original.external_metabolite_indices()]
#print(external_met_IDs_original)
#external_model_sids_original = [species.id for species in intermediate_cmod.species if species.compartment == 'e']
#print(external_model_sids_original)

#print(set(external_met_IDs_original) - set(external_model_sids_original))
#print(set(external_model_sids_original) - set(external_met_IDs_original))

#print(set(external_met_IDs) - set(external_met_IDs_original))
#print(set(external_met_IDs_original) - set(external_met_IDs))

"""Calculate ECMs and/or EFMs"""
# For genome-scale models we cannot yet calculate ECMs, then we should calculate them for only the active subnetwork,
# or only for the model with external metabolites that are ignored.
# TODO: Make calculation stop when it takes too long, and give error message
# For now I manually remove some models defined at top of the file
#DROP_MODELS = [1,2,3]
list_model_dicts_remember = list_model_dicts.copy()
list_model_dicts = [model_dict for index, model_dict in enumerate(list_model_dicts) if index not in DROP_MODELS]

# Calculate ECMs for all models, and add it to the model dictionaries
for model_dict in list_model_dicts:
    ecms_df, full_network_ecm, efms_df, ecms_from_efms = get_ecm_df(result_dir, model_dict['model_path'], infos_obj,
                                                                    hide_indices=model_dict['hide_metabolites'],
                                                                    to_be_tagged=reactions_to_tag,
                                                                    print_results=True,
                                                                    get_EFMs=model_dict['calc_efms'],
                                                                    external_compartment=EXTERNAL_COMPARTMENT)
    model_dict.update(
        {'ecms_df': ecms_df, 'network': full_network_ecm, 'efms_df': efms_df, 'ecms_from_efms': ecms_from_efms})

# Make a list with for each model how many ECMs were found
list_n_ECMs = [(model_dict['ecms_df'].shape[1], len(model_dict['model'].getSpeciesIds()),
                len(model_dict['model'].getReactionIds())) for model_dict in list_model_dicts]
# Print the number of metabolites, reactions and ECMs
for model_ind, counted_model in enumerate(list_n_ECMs):
    print('Model %d has %d metabolites, %d reactions and %d ECMs' % (
        model_ind, counted_model[1], counted_model[2], counted_model[0]))

# TODO: Add annotation to Matlab-file
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

"""Save the data in list_model_dicts. Also save important parameters of this session. Especially when needed for plotting."""
# Save list_model_dicts to file. pickle? To be able to replot/revisit the calculations.

# Store data
# with current data_time in name ? then it cannot be overwritten
import datetime
# TODO: Remove these print-statements or add annotation
now = datetime.datetime.now()
print("Current date and time : ")
print(now.strftime("%Y-%m-%d %H:%M:%S"))
print(now.strftime("%Y%m%d_%H%M%S"))

filename = 'list_model_dicts' + now.strftime("_%Y%m%d_%H%M%S") + '.pickle'
with open(os.path.join(result_dir, filename), 'wb') as handle:
    pickle.dump(list_model_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Load data
# Check specific date_time !!
# filename = 'list_model_dicts_20201006_145107.pickle'
# with open(os.path.join(result_dir, '..', '..', 'e_coli_core_2', 'EFM_yield_analysis', filename), 'rb') as handle:
#    list_model_dicts = pickle.load(handle)

"""Plot these costs"""
# We need the objective value to scale the axes of our plots
objective_val = intermediate_cmod.getObjFuncValue()
# We make 2-D cost plots only (3-D is too confusing). Since each constrained metabolite brings costs, each metabolite
# needs to be plotted on some axes. We thus make ceil(n_cons/2) figures
constrained_ext_metabs = [cons_dict['ext_metab'] for cons_dict in list_model_dicts[0]['infos_cons']]
for ind in range(0, len(constrained_ext_metabs), 2):
    # print(ind)
    # Select constraints that will be plotted in this figure
    cons_for_fig = constrained_ext_metabs[ind:max(ind + 2, len(constrained_ext_metabs))]
    # Plot the results of the different models (with different (sub)networks) in one plot
    plot_different_approaches_one_figure(list_model_dicts, infos_obj,  # list_model_dicts_remember[0]['infos_obj'],
                                         infos_cons,  # list_model_dicts_remember[0]['infos_cons'],
                                         obj_val=objective_val,
                                         result_dir=result_dir, cons_IDs=cons_for_fig)

# Now make for each approach for which we calculated the activities of the ECMs a different plot in which the usage of
# the ECMs is shown with vectors.
if len(constrained_ext_metabs) > 1:
    cons_ID_combi_list = get_cons_ID_combinations(constrained_ext_metabs)
else:
    cons_ID_combi_list = None

if cons_ID_combi_list:
    for model_dict in list_model_dicts:
        if model_dict['get_activities']:
            print('Plotting the cost vectors including usage for model %s' % model_dict['model_name'])
            for cons_ID_combi in cons_ID_combi_list:
                plot_costs(model_dict, infos_obj, infos_cons,  # model_dict['infos_obj'], model_dict['infos_cons'],
                           cons_IDs=cons_ID_combi, obj_val=objective_val,
                           show_active=True, result_dir=result_dir)
else:
    for model_dict in list_model_dicts:
        if model_dict['get_activities']:
            print('Plotting the cost vectors including usage for model %s' % model_dict['model_name'])
            plot_costs(model_dict, infos_obj, infos_cons,  # model_dict['infos_obj'], model_dict['infos_cons'],
                       cons_IDs=constrained_ext_metabs, obj_val=objective_val,
                       show_active=True, result_dir=result_dir)
#
# for model_dict in list_model_dicts:
#     if model_dict['get_activities']:
#         print('Plotting the cost vectors including usage for model %s' % model_dict['model_name'])
#         plot_costs(model_dict, infos_obj, infos_cons, #model_dict['infos_obj'], model_dict['infos_cons'],
#                    cons_IDs=constrained_ext_metabs, obj_val=objective_val,
#                    show_active=True, result_dir=result_dir)

# TODO: Make a function that finds an EFM corresponding to an ECM
# TODO: Apply this function to the resulting ECMs for the hidden metabolites, to get a minimal network needed for these constraints
# TODO: Apply this function (if feasible) for all ECMs, to get a minimal FBA-network

# Find some EFM that corresponds to each of the active ECMs in the hidden-network
# for model_dict in list_model_dicts:
    if model_dict['get_relevant_efms'] and model_dict['get_activities']:
        relevant_efms_df, full_relevant_ecms_df = find_associated_efms(model_dict['model'], model_dict['table_cons_df'],
                                                                       model_dict['ecms_df'],
                                                                       infos_obj + infos_cons, model_dict['model_path'],
                                                                       external_compartment=EXTERNAL_COMPARTMENT)
        relevant_efms_df.to_csv(
            os.path.join(
                result_dir, "efms_corresponding_to_hide_ecms" + model_dict['model_name'] + ".csv"), index=False)
        full_relevant_ecms_df.to_csv(
            os.path.join(
                result_dir, "full_ecms_corresponding_to_hide_ecms" + model_dict['model_name'] + ".csv"), index=False)

""" problem solving """
# print active ECMs
# for model_dict in list_model_dicts:
#     print('Model =', model_dict['model_name'])
#     print('Activity per ECM for objective(s)')
#     print(model_dict['table_obj_df'][model_dict['table_obj_df']['active']!=0.])
#     print('Activity per ECM for the active constraint(s)')
#     print(model_dict['table_cons_df'][model_dict['table_cons_df']['active']!=0.])
#     print('Glucose flux:', sum(model_dict['table_cons_df']['active']*model_dict['table_cons_df']['M_glc__D_e']))
#     print('Oxygen flux:', sum(model_dict['table_cons_df']['active']*model_dict['table_cons_df']['M_o2_e']))
#     print()
#     print('Active ECMs and used metabolites')
#     print(get_active_ecms(model_dict))
#     print()