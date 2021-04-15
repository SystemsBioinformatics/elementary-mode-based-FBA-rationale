import itertools

from helpers_EFM_FBA import *


"""
NOTE: Run this script three times with different uptake bounds for oxygen uptake (-15, -10, -5) to generate three cost 
vector plots that are the base for fig 3.
"""


"""CONSTANTS"""
model_name = "iJR904"

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

################ Adapt oxygen uptake reaction here #####################
oxygen_reaction.setLowerBound(-10.)  # -15, -10, -5-> give two active constraints; -1000 -> give 1 active constraint
########################################################################

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
cbm.doFBAMinSum(intermediate_cmod)

"""Create list with a dictionary in which all info will be stored for each model"""
list_model_dicts = []

# Create dictionary for original model
list_model_dicts.append(
    {'model': intermediate_cmod, 'model_name': 'original_network', 'calc_efms': False, 'get_activities': True,
     'get_relevant_efms': False, 'drop_tag': 'full'})

if ('active' not in DROP_MODEL_TAGS) or ('active_hidden' not in DROP_MODEL_TAGS):
    # Deletes all non active reactions (determined by seeing if the flux is lower than a tolerance)
    sub_model_name = 'active_subnetwork'
    active_model_path = os.path.join(sbml_dir, sub_model_name + ".xml")
    if os.path.exists(active_model_path) and LOAD_IF_AVAILABLE:  # Checks if this model was already made
        try:
            intermediate_cmod_active = cbm.readSBML3FBC(active_model_path)
        except:  # use version 2
            intermediate_cmod_active = cbm.readSBML2FBA(active_model_path)
    else:
        intermediate_cmod_active = delete_non_active_network(intermediate_cmod, which_zeros='flux_zero', zero_tol=1e-15,
                                                             opt_tol=1e-8)  # 15, 8

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

if 'hidden' not in DROP_MODEL_TAGS:
    """Get relevant metabolites for the cost-calculations, find indices of external metabs that can be ignored"""
    # We only have to focus on conversions of metabolites that are somehow active in a constraint
    relevant_metabs = [info_dict['ext_metab'] for info_dict in infos_obj + infos_cons]

    # The following finds all indices of external metabolites that can be ignored in the conversion analysis
    hide_indices = find_hide_indices(full_model_path, to_be_tagged=reactions_to_tag, focus_metabs=relevant_metabs,
                                     use_external_compartment=USE_EXTERNAL_COMPARTMENT)
# By default no metabolites are hidden
for model_dict in list_model_dicts:
    model_dict['hide_metabolites'] = []

if 'hidden' not in DROP_MODEL_TAGS:
    """Create another dictionary for a model in which some metabolites are hidden"""
    list_model_dicts.append(
        {'model': intermediate_cmod, 'model_name': 'original_with_hidden_metabolites', 'calc_efms': False,
         'get_activities': True, 'hide_metabolites': hide_indices, 'get_relevant_efms': True,
         'model_path': os.path.join(sbml_dir, "original_network.xml"), 'drop_tag': 'hidden'})

if 'active_hidden' not in DROP_MODEL_TAGS:
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

# We remove some models defined at top of the file
list_model_dicts_remember = list_model_dicts.copy()
list_model_dicts = [model_dict for model_dict in list_model_dicts if model_dict['drop_tag'] not in DROP_MODEL_TAGS]

# Load ECMs file E.coli iJR
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
if len(constrained_ext_metabs) > 1:
    if model_name == "e_coli_core":
        # all possible combinations of constraints
        cons_ID_combi_list = list(itertools.combinations(constrained_ext_metabs, 2))
    else:
        cons_ID_combi_list = get_cons_ID_combinations(constrained_ext_metabs, list_model_dicts[0])
else:
    cons_ID_combi_list = None

# Now make for each approach for which we calculated the activities of the ECMs a different plot in which the usage of
# the ECMs is shown with vectors.
if cons_ID_combi_list:
    for model_dict in list_model_dicts:
        if model_dict['get_activities']:
            print('Plotting the cost vectors including usage for model %s' % model_dict['model_name'])
            for cons_ID_combi in cons_ID_combi_list:
                print(cons_ID_combi)
                plot_costs_flux(model_dict, infos_obj, infos_cons,  # model_dict['infos_obj'], model_dict['infos_cons'],
                           cons_IDs=cons_ID_combi, obj_val=objective_val,
                           show_active=True, result_dir=result_dir, squared_plot=False)

"""Find corresponding EFM(s) for active ECM(s)"""
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
            printable_ecms = full_relevant_ecms_df.sort_values('M_glc__D_e')/2
            print_ecms_direct(np.transpose(printable_ecms.values), full_relevant_ecms_df.columns)


""" Supplemental figure A4 on ECM product formation per oxygen uptake """
# Select relevant ECMs (=biomass producing ECMs) and relevant metabolites (metabolites for which stoichiometry is not the same in all relevant ECMs)
relevant_metabolites = ["M_ac_e", "M_co2_e", "M_etoh_e", "M_for_e", "M_glc__D_e", "M_glyald_e", "M_glyclt_e", "M_h2o_e", "M_h_e", "M_o2_e", "M_succ_e", "objective"]
relevant_ECMs = full_relevant_ecms_df[full_relevant_ecms_df['objective']!=0.]
sorted_ecms = relevant_ECMs.sort_values('M_o2_e')
sorted_ecms[relevant_metabolites]
sorted_ecms['ox_per_glc'] = sorted_ecms['M_o2_e']/sorted_ecms['M_glc__D_e']
sorted_ecms['Y_biomass_glc'] = sorted_ecms['objective']/-sorted_ecms['M_glc__D_e']

# Select metabolites to be plotted
plot_metabolites = ["M_ac_e", "M_co2_e", "M_etoh_e", "M_for_e", "M_glyclt_e", "M_succ_e"]

# Make plot: product_formation_per_ox_per_glc
fig, ax = plt.subplots(1, 1) #, figsize=(7, 4))
for metabolite in plot_metabolites:
    # plot rates per 0.5 unit biomass
    plt.plot(sorted_ecms['ox_per_glc'], sorted_ecms[metabolite]/2,
                label=intermediate_cmod.getSpecies(metabolite).getName().split(" ")[0],
                marker='.',
                linestyle='-')
plt.legend(loc='upper right', ncol=2)
plt.ylabel('Product formation per 0.5 unit biomass')
plt.xlabel('Oxygen per glucose')
plt.xticks(rotation = 90)
plt.tight_layout()
plt.savefig(os.path.join(result_dir, "product_formation_per_ox_per_glc" + ".png"))
plt.savefig(os.path.join(result_dir, "product_formation_per_ox_per_glc" + ".pdf"))
plt.savefig(os.path.join(result_dir, "product_formation_per_ox_per_glc" + ".svg"))
plt.show()