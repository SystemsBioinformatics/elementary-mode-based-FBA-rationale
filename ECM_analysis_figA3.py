import itertools

from helpers_EFM_FBA import *

"""CONSTANTS"""
model_name = "Lactobacillus_plantarum_WCFS1_Official_23_May_2019_18_45_01"

# Define directories for storing: 1) newly made models, 2) results such as figures, and 3) for finding models
sbml_dir = os.path.join(os.getcwd(), "data", model_name + "_2", 'models', 'sbml')
result_dir = os.path.join(os.getcwd(), "data", model_name + "_2", "EFM_yield_analysis")
model_dir = os.path.join("models")
LOAD_IF_AVAILABLE = False  # If this is true, newly created models will be loaded if the program runs a second time
N_KOS = 0  # number of knockout models that will be created

# The complete model is sometimes to large to calculate ECMs. This is therefore dropped manually
# adjust DROP_MODELS based on this.
DROP_MODEL_TAGS = ['full', 'active', 'fva', 'hidden', 'ko']  # ['full','active','fva','hidden','active_hidden' ,'ko']
USE_EXTERNAL_COMPARTMENT = None
ONLY_RAYS = False # True or False
EXCECPTIONS = [] # reaction IDs that shouldn't be deleted in active model
N_COLS = 3 # number of columns in plot_costs_ECMS_one_figure
SQUARED = True # forced squared plot for cost vector plot

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
if model_name == "Lactobacillus_plantarum_WCFS1_Official_23_May_2019_18_45_01":
    EXCECPTIONS = ['R_PYDXK2'] # no deletion of reaction in active model
    PLANTARUM_BOUNDS = 'teusink2009_medium' # Choose from 'default', 'teusink2009', 'teusink2009_medium'
    if PLANTARUM_BOUNDS == 'default':
        for rid in intermediate_cmod.getReactionIds():
            print(intermediate_cmod.getReactionBounds(rid))
            intermediate_cmod.setReactionBounds(rid, intermediate_cmod.getReactionLowerBound(rid) * 100.,
                                                     intermediate_cmod.getReactionUpperBound(rid) * 100.)
    elif PLANTARUM_BOUNDS == 'teusink2009':
        filepath = os.path.join(result_dir, "..", "Teusink2009bounds.xlsx")
        teusink2009_bounds = pd.read_excel(filepath, index_col=1)

        # Set biomass objective
        intermediate_cmod.setReactionBounds('R_biomass_LPL_RETB_t576_NoATP', 0., 0.)
        intermediate_cmod.setReactionBounds('R_biomass_LPL60', 0., 999999.)
        intermediate_cmod.createObjectiveFunction('R_biomass_LPL60')

        # make reaction reversible because of infeasible results & should be reversible transporter
        reaction = intermediate_cmod.getReaction("R_THRt2r")
        reaction.reversible=True
        reaction.setLowerBound(-999999.)

        # Set bounds exchange reactions
        for rid in intermediate_cmod.getReactionIds("R_EX_"):
            intermediate_cmod.setReactionBounds(rid, 0., 999999.)
        for rid in teusink2009_bounds['RID']:
            intermediate_cmod.setReactionLowerBound(rid, teusink2009_bounds[teusink2009_bounds['RID']==rid]['LOWER BOUND'][0])
            intermediate_cmod.setReactionUpperBound(rid, teusink2009_bounds[teusink2009_bounds['RID']==rid]['UPPER BOUND'][0])

        for rid in teusink2009_bounds['RID']:
            if teusink2009_bounds[teusink2009_bounds['RID']==rid]['LOWER BOUND'][0] > 0:
                intermediate_cmod.setReactionLowerBound(rid,0.)
                intermediate_cmod.setReactionUpperBound(rid,999999.)
            elif teusink2009_bounds[teusink2009_bounds['RID']==rid]['LOWER BOUND'][0] == 0:
                intermediate_cmod.setReactionLowerBound(rid,0.)
                intermediate_cmod.setReactionUpperBound(rid,999999.)
            else:# LB < 0.
                #intermediate_cmod.setReactionLowerBound(rid,-999999.)
                if teusink2009_bounds[teusink2009_bounds['RID']==rid]['UPPER BOUND'][0] <= 0:
                    intermediate_cmod.setReactionUpperBound(rid,0.)
                else: #UB > 0.
                    intermediate_cmod.setReactionUpperBound(rid,999999.)

        # Test if higher ATP maintenance requirement leads to same growth rate as in paper Teusink 2009. Yes. 0.26/h
        #intermediate_cmod.setReactionBounds("R_ATPM", 3.94, 3.94)
        intermediate_cmod.setReactionBounds("R_ATPM", 0.29, 0.29)

        # Robustness_analysis_changes
        newteusink2009_bounds = pd.read_excel(filepath, sheet_name='Robustness_analysis_changes', index_col=1)
        for rid in newteusink2009_bounds['RID']:
            intermediate_cmod.setReactionUpperBound(rid, newteusink2009_bounds[newteusink2009_bounds['RID']==rid]['UPPER BOUND'][0])
            intermediate_cmod.setReactionLowerBound(rid, newteusink2009_bounds[newteusink2009_bounds['RID']==rid]['LOWER BOUND'][0])

        intermediate_cmod.setReactionBounds("R_EX_o2_e",
                                            teusink2009_bounds[teusink2009_bounds['RID'] == "R_EX_o2_e"]['LOWER BOUND'][0],
                                            teusink2009_bounds[teusink2009_bounds['RID'] == "R_EX_o2_e"]['UPPER BOUND'][0])

        # Change all UB 999999 to inf and all LB -999999 to -inf
        reactions = cmod.getReactionIds()
        for rid in reactions:
            if intermediate_cmod.getReactionLowerBound(rid) == -999999:
                intermediate_cmod.setReactionLowerBound(rid, -np.inf)
            if intermediate_cmod.getReactionUpperBound(rid) == 999999.:
                intermediate_cmod.setReactionUpperBound(rid, np.inf)

    elif PLANTARUM_BOUNDS == 'teusink2009_medium':
        filepath = os.path.join(result_dir, "..", "Teusink2009bounds.xlsx")
        teusink2009_bounds = pd.read_excel(filepath, index_col=1, sheet_name='Medium')

        # Set biomass objective
        intermediate_cmod.setReactionBounds('R_biomass_LPL_RETB_t576_NoATP', 0., 0.)
        intermediate_cmod.setReactionBounds('R_biomass_LPL60', 0., 999999.)
        intermediate_cmod.createObjectiveFunction('R_biomass_LPL60')

        # make reaction reversible because should be a reversible transporter
        reaction = intermediate_cmod.getReaction("R_THRt2r")
        reaction.reversible = True
        reaction.setLowerBound(-999999.)

        # Set bounds exchange reactions
        for rid in intermediate_cmod.getReactionIds("R_EX_"):
            intermediate_cmod.setReactionBounds(rid, 0., 999999.)
        for rid in teusink2009_bounds['RID']:
            intermediate_cmod.setReactionLowerBound(rid,
                                                    teusink2009_bounds[teusink2009_bounds['RID'] == rid]['LOWER BOUND'][
                                                        0])
            intermediate_cmod.setReactionUpperBound(rid,
                                                    teusink2009_bounds[teusink2009_bounds['RID'] == rid]['UPPER BOUND'][
                                                        0])

        # Robustness_analysis_changes
        newteusink2009_bounds = pd.read_excel(filepath, sheet_name='Robustness_analysis_changes', index_col=1)
        for rid in newteusink2009_bounds['RID']:
            intermediate_cmod.setReactionUpperBound(rid, newteusink2009_bounds[newteusink2009_bounds['RID'] == rid][
                'UPPER BOUND'][0])
            intermediate_cmod.setReactionLowerBound(rid, newteusink2009_bounds[newteusink2009_bounds['RID'] == rid][
                'LOWER BOUND'][0])

        intermediate_cmod.setReactionBounds("R_EX_o2_e",
                                            teusink2009_bounds[teusink2009_bounds['RID'] == "R_EX_o2_e"]['LOWER BOUND'][
                                                0],
                                            teusink2009_bounds[teusink2009_bounds['RID'] == "R_EX_o2_e"]['UPPER BOUND'][
                                                0])

        # Change all UB 999999 to inf and all LB -999999 to -inf
        reactions = cmod.getReactionIds()
        for rid in reactions:
            if intermediate_cmod.getReactionLowerBound(rid) == -999999:
                intermediate_cmod.setReactionLowerBound(rid, -np.inf)
            if intermediate_cmod.getReactionUpperBound(rid) == 999999.:
                intermediate_cmod.setReactionUpperBound(rid, np.inf)

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
    except:  # use version 2
        intermediate_cmod_active = cbm.readSBML2FBA(active_model_path)
else:
    intermediate_cmod_active = delete_non_active_network(intermediate_cmod, which_zeros='flux_zero', zero_tol=1e-20,
                                                         opt_tol=1e-8, exceptions=EXCECPTIONS) # 15, 8

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

"""Get relevant metabolites for the cost-calculations, find indices of external metabs that can be ignored"""
# We only have to focus on conversions of metabolites that are somehow active in a constraint
relevant_metabs = [info_dict['ext_metab'] for info_dict in infos_obj + infos_cons]
# The following finds all indices of external metabolites that can be ignored in the conversion analysis
hide_indices = find_hide_indices(full_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs,
                                 use_external_compartment=USE_EXTERNAL_COMPARTMENT)
# By default no metabolites are hidden
for model_dict in list_model_dicts:
    model_dict['hide_metabolites'] = []

"""Create another dictionary for a model in which some metabolites are hidden"""
list_model_dicts.append(
    {'model': intermediate_cmod, 'model_name': 'original_with_hidden_metabolites', 'calc_efms': False,
     'get_activities': True, 'hide_metabolites': hide_indices, 'get_relevant_efms': True,
     'model_path': os.path.join(sbml_dir, "original_network.xml"), 'drop_tag': 'hidden'})

"""Create another dictionary for the active model with some metabolites hidden"""
active_hide_indices = find_hide_indices(active_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs,
                                        use_external_compartment=USE_EXTERNAL_COMPARTMENT)
list_model_dicts.append(
    {'model': intermediate_cmod_active, 'model_name': 'active_network_with_hidden_metabolites', 'calc_efms': False, #False
     'get_activities': True, 'hide_metabolites': active_hide_indices, 'get_relevant_efms': True,
     'model_path': active_model_path,  'drop_tag': 'active_hidden'})

"""Calculate ECMs and/or EFMs"""
# For genome-scale models we cannot yet calculate ECMs, then we should calculate them for only the active subnetwork,
# or only for the model with external metabolites that are ignored.
# TODO: Make calculation stop when it takes too long, and give error message
# For now I manually remove some models defined at top of the file
list_model_dicts_remember = list_model_dicts.copy()
list_model_dicts = [model_dict for model_dict in list_model_dicts if model_dict['drop_tag'] not in DROP_MODEL_TAGS]

# Calculate ECMs for all models, and add it to the model dictionaries
for model_dict in list_model_dicts:
    ecms_df, full_network_ecm, efms_df, ecms_from_efms = get_ecm_df(result_dir, model_dict['model_path'], infos_obj,
                                                                    hide_indices=model_dict['hide_metabolites'],
                                                                    to_be_tagged=reactions_to_tag,
                                                                    print_results=False, only_rays=ONLY_RAYS,
                                                                    get_EFMs=model_dict['calc_efms'],
                                                                    use_external_compartment=USE_EXTERNAL_COMPARTMENT)
    model_dict.update(
        {'ecms_df': ecms_df, 'network': full_network_ecm, 'efms_df': efms_df, 'ecms_from_efms': ecms_from_efms})

# Make a list with for each model how many ECMs were found
list_n_ECMs = [(model_dict['ecms_df'].shape[1], len(model_dict['model'].getSpeciesIds()),
                len(model_dict['model'].getReactionIds())) for model_dict in list_model_dicts]
# Print the number of metabolites, reactions and ECMs
for model_ind, counted_model in enumerate(list_n_ECMs):
    print('Model %d has %d metabolites, %d reactions and %d ECMs' % (
        model_ind, counted_model[1], counted_model[2], counted_model[0]))

"""Select only the relevant part of the ECM information. So, only the consumption/production flux of metabolites that
    are associated to a flux bound."""
for model_dict in list_model_dicts:
    model_dict['infos_cons'] = infos_cons
    model_dict['infos_obj'] = infos_obj
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

if 'Lactobacillus_plantarum_WCFS1_Official_23_May_2019_18_45_01' in model_name:
    # Get cons_IDs in combinations in order of abs(scaled reduced cost)
    cons_ID_combi_list = get_cons_ID_combinations(constrained_ext_metabs, model_dict)
    # In latter case: adjust vector_focus to vector_focus = False in plot_costs_ECMS_one_figure()

# Now make for each approach for which we calculated the activities of the ECMs a different plot in which the usage of
# the ECMs is shown with vectors.
if cons_ID_combi_list:
    for model_dict in list_model_dicts:
        if model_dict['get_activities']:
            print('Plotting the cost vectors including usage for model %s' % model_dict['model_name'])

            # Plot costs of ECMS in constraint space plot
            plot_costs_ECMS_one_figure(model_dict, infos_cons, cons_ID_combi_list=cons_ID_combi_list,
                                       result_dir=result_dir, n_cols=N_COLS, vector_focus=False)
