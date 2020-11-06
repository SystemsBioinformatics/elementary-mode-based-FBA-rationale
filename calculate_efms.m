function efm_struct = calculate_efms(efmtool_folder,model_name,xml_file)
%initCobraToolbox
%model_name = 'toy_model';
%model_name = 'core_memesa_model.l3fbc';
%model_name = 'iTM686.light';
%model_name = 'iAF1260';
%model_name = 'iYO844';
%model_name = 'iMM904';
%model_name = 'NCDO712efmtool2';
%model_name = 'e_coli_core_2';
%model_name = 'iAF1260_2';
%model_name = 'iTM686.light_2';
%model_name = 'iMM904_2';
%model_name = 'iYO844_2';
%model_name = 'core_memesa_model.l3fbc_2';
%folder_name = sprintf(model_name)

%fileID = fopen('folder_names.txt');
%foldernames=textscan(fileID,'%s');
%fclose(fileID);

%efmtool_folder=foldernames{1}{1};
%model_name=foldernames{1}{2};

% go to directory
dirname=cd;
cd(sprintf('data'));
cd(model_name);
cd(sprintf('EFM_yield_analysis'));

% load model
%model = readCbModel('lumped_model.xml');
model = readCbModel(xml_file); %, 'fileType','SBML','defaultBound', 999999); model = readSBML('reduced_model.xml');
% 
model.S = full(model.S);

model.rev = [];
for K = 1:size(model.lb)
    if model.lb(K) < 0
        model.rev = [model.rev; 1];
    else
        model.rev = [model.rev; 0];
    end
end
%% 

cd(efmtool_folder);

%rxnname = ['R26']; %804   'R_BiomassAuto'
%opts = CreateFluxModeOpts('enforce', 'R25')

%%
opts = CreateFluxModeOpts('arithmetic', 'fractional');
mnet = CalculateFluxModes(model.S, model.rev, model.mets, model.rxns,opts); %, opts);
%           - stoich           the stoichiometric matrix
%           - reversibilities  reaction reversibilities
%                              0/1 for irreversible/reversible
%           - mnames           metabolite names (cell array)
%           - rnames           reaction names (cell array)

%%

cd(dirname);
cd(sprintf('data'));
cd(sprintf(model_name));
cd(sprintf('EFM_yield_analysis'));

%%

efm_struct.efms = mnet.efms;
efm_struct.reacNames=mnet.reactionNames;
%reacNames = mnet.reactionNames;
%T = table(reacNames', efms);
%filename = 'efmdata.xlsx';
%if exist(filename, 'file') == 2
%    delete(filename)
%end
%writetable(T,filename)



%fid = fopen('test_file','w');
%for k=1:length(reacNames)
%   fprintf(fid,'%s\t%.2f\t%.2f\t%.2f\n',reacNames{k},efms(k,:));
%end
%fclose(fid);
