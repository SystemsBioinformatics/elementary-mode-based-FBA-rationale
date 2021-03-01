function efms = calculate_efms_wo_SBML(efmtool_folder, result_dir)

cd(result_dir)

% Convert stoich_list to matrix
stoich_matrix = readtable('stoich_matrix.csv', 'ReadVariableNames',false);
reversibilities = readtable('reversibilities.csv', 'ReadVariableNames',false);

stoich_matrix = table2array(stoich_matrix);
reversibilities = table2array(reversibilities);

cd(efmtool_folder);

opts = CreateFluxModeOpts('arithmetic', 'double');
mnet = CalculateFluxModes(stoich_matrix, reversibilities,opts);

efms =  mnet.efms;