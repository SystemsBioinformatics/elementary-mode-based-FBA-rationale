from helpers_EFM_FBA import *

model_dir = os.path.join("models")
model_list = ["iJR904",
              "Lactobacillus_plantarum_WCFS1_Official_23_May_2019_18_45_01",
              "MG1363_20190628",
              "iAF1260b",
              "iMM904",
              "iNJ661"]
for model_name in model_list:
    # load model
    try:
        cmod = cbm.readSBML3FBC(os.path.join(model_dir, model_name + ".xml"))
        # cmod.setNotes(cmod.getNotes().encode(encoding='ascii', errors='ignore'))
    except:  # use version 2
        cmod = cbm.readSBML2FBA(os.path.join(model_dir, model_name + ".xml"))

    # do FBA with scaled reduced cost
    cbm.analyzeModel(cmod, with_reduced_costs='scaled')

    # find scaled reduced costs
    nzrc_dictionaries_scaled, n_objectives_scaled = get_nzrc(cmod)

    # report output
    print('Number of nzrc reactions for model:')
    print(model_name)
    print(len(nzrc_dictionaries_scaled))