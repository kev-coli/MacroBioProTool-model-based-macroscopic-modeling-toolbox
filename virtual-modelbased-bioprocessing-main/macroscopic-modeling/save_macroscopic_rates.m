addpath('C:\Program Files\Mosek\10.2\toolbox\r2017aom')
name_save_file = "Test2";
directory_file = strcat("../models/",name_save_file,"/");
pathway_matlab_file_model = strcat(directory_file,name_save_file,".mat");
name_excel_saving_file = strcat(directory_file,strcat(name_save_file,".xls"));
load(pathway_matlab_file_model);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the EFM column generation again with fixed subset of EFMs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qext_train = data_qext.qext_train';
qext_predict = data_qext.qext_predict';
nb_data_per_training_media = data_qext.nb_data_per_training_media;
nb_data_per_prediction_media = data_qext.nb_data_per_prediction_media;
options_qp = optimoptions('quadprog','Display','none'); % 
average_data_per_condition = options_data.average_data_per_condition;

%% Load the necessary data
index_training_media = data_qext.index_training_media;
index_prediction_media = data_qext.index_prediction_media;
use_mosek = options_EFM.use_mosek;
ic = data_qext.ic;
qext_train = data_qext.qext_train';
qext_predict = data_qext.qext_predict';
Ameas = data_stoich.Ameas; 
Aint = data_stoich.Aint;
Dirrev = data_stoich.Dirrev; 
n_v = data_stoich.n_v;
n_v_irrev = data_stoich.n_v_irrev;
n_mets_int = data_stoich.n_mets_int;
[n_qext,N_train] = size(qext_train);
[~,N_predict] = size(qext_predict);

tolerance_EFM = options_EFM.tolerance_EFM;
computation_method_EFM = options_EFM.computation_method_EFM;
reduction_of_EFM = options_EFM.reduction_of_EFM;
factor_error_reduction = options_EFM.factor_error_reduction;
normalization_matrix = data_qext.normalization_matrix;
Ameas = normalization_matrix*Ameas;

Aint_sp = sparse(Aint); % Convert to sparse matrix in order to increase accuracy for optimization
Ameas_sp = sparse(Ameas); % Convert to sparse matrix in order to increase accuracy for optimization
if(use_mosek == 0)
  H = Ameas'*Ameas + 1e-10*eye(n_v); % Quadratic matrix for the Matlab function quadprog
  H = (H+H')/2; % This line guarantees that the quadratic matrix is symmetric
else
  H = [];  
end
Aeq_qp = Aint_sp; % Matrix used to implement the equality constraints Aint*v = 0
beq_qp = sparse(n_mets_int,1); % Vector used to implement the equality constraints Aint*v = 0
Aineq_qp = -sparse(Dirrev);  % Matrix used to implement the inequality constraints Dirrev*v >= 0 
bineq_qp = sparse(n_v_irrev,1); % Vector used to implement the equality constraints Dirrev*v >= 0 
qext_proj = zeros(n_qext,N_train); % Matrix containing the projected cell specific rate data onto the realisable set
for kk = 1:N_train  
  Yk = qext_train(:,kk); fk = sparse(-Ameas_sp'*Yk);
  if(use_mosek) % solve the quadratic problem with MOSEK
    v_opt_kk = quadprog_mosek_as_second_order_programming(Ameas_sp,Yk,Aineq_qp,bineq_qp,Aeq_qp,beq_qp,-Inf*ones(n_v,1),Inf*ones(n_v,1));  
  else  % solve the quadratic problem with quadprog
    v_opt_kk = quadprog(H,fk,Aineq_qp,bineq_qp,Aeq_qp,beq_qp,[],[],[],options_qp);
  end
  qext_proj(:,kk) = Ameas_sp*v_opt_kk; 
end
Y = qext_proj; K = Ameas; 
[~,~,~,Macro_rates_train] = master_and_subproblem(qext_proj,K,EFMs,Aint,Dirrev,use_mosek);
if(N_predict > 0)
  [~,~,~,Macro_rates_predict] = master_and_subproblem(qext_predict,K,EFMs,Aint,Dirrev,use_mosek);
else
  Macro_rates_predict = [];    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Reconstruct kinetic models of the macroscopic rates       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mets_not_included_in_Monod_kinetics = options_kinetic_identification.mets_not_included_in_Monod_kinetics;
mets_meas_in_cext = data_cext.mets_meas_in_cext;
cext = data_cext.cext; % load the concentration data
[~,n_cext] = size(cext);
[n_macro_rates,~] = size(Macro_rates_train); 
cext_train = data_cext.cext_train; % load the concentration TRAINING data
cext_predict = data_cext.cext_predict; % load the concentration PREDICTION data
cext_train(cext_train < 1e-4) = 1e-4; % Replace zero concentration by a small value. This prevents numerical instability during modeling
cext_predict(cext_predict < 1e-4) = 1e-4; % Replace zero concentration by a small value. This prevents numerical instability during modeling
input = cext_train;
input_predict = cext_predict;
  
% Exclude the concentration data which are not included in the kinetics (decided by the user with the cell vector mets_not_included_in_Monod_kinetics,mets_meas_in_cext)
[~,locc] = ismember(mets_not_included_in_Monod_kinetics,mets_meas_in_cext);
locc = nonzeros(locc); ind_cext_in_kinetic_identification = setdiff(1:n_cext,locc);
input_train = input(:,ind_cext_in_kinetic_identification); % input vector for the Monod model. Each column corresponds to a concentration of a metabolite which possibly participate in the macroscopic reactions
if(N_predict > 0)
  input_predict = input_predict(:,ind_cext_in_kinetic_identification);
end
n_inputs = size(input,2);


activation_parameters = Kinetic_parameters.activation_parameters;
inhibition_parameters = Kinetic_parameters.inhibition_parameters;
max_reaction_rate = Kinetic_parameters.max_reaction_rate;
N_train = size(input,1); N_predict = size(input_predict,1);
w_model_train = ones(n_macro_rates,N_train); w_model_predict = ones(n_macro_rates,N_predict);
for j = 1:n_macro_rates
   act_j = activation_parameters(j,:)'; inh_j = inhibition_parameters(j,:)'; alp_j = max_reaction_rate(j,1);
   for k = 1:n_inputs
     w_model_train(j,:) = w_model_train(j,:).*(input(:,k)')./(input(:,k)' + act_j(k))./(1 + inh_j(k)*input(:,k)');
     if(N_predict > 0)
       w_model_predict(j,:) = w_model_predict(j,:).*(input_predict(:,k)')./(input_predict(:,k)' + act_j(k))./(1 + inh_j(k)*input_predict(:,k)');
     end
   end
   w_model_train(j,:) = w_model_train(j,:)*alp_j;
   if(N_predict > 0)
     w_model_predict(j,:) = w_model_predict(j,:)*alp_j;
   end
end

%% COMPLETING THE EXCEL FILES

media = data_qext.media;

% we save the cell specific rate computed from the kinetic model
  index_training_media = data_qext.index_training_media; 
  index_prediction_media = data_qext.index_prediction_media; 
  N = data_qext.N;
  nb_data_per_training_media = data_qext.nb_data_per_training_media;
  nb_data_per_prediction_media = data_qext.nb_data_per_prediction_media;
  ic = data_qext.ic;
  days_column = data_qext.days_column;
  cell_sheet9 = cell(3+n_macro_rates,N);
  kl_index = 1;
  for k = 1:length(nb_data_per_training_media)
    nb_data_in_training_media_k = nb_data_per_training_media(k);  
    kk = index_training_media(k);
    ind_days_media_k = (ic == kk);
    all_days_for_media_k = days_column(ind_days_media_k);
    for j = 1:nb_data_in_training_media_k
      cell_sheet9(1,kl_index) = {"Training"};
      cell_sheet9(2,kl_index) = media(index_training_media(k));
      cell_sheet9(3,kl_index) = all_days_for_media_k(j);
      kl_index =  kl_index + 1;
    end
  end
  for k = 1:length(nb_data_per_prediction_media)
    nb_data_in_prediction_media_k = nb_data_per_prediction_media(k);  
    kk = index_prediction_media(k);
    ind_days_media_k = (ic == kk);
    all_days_for_media_k = days_column(ind_days_media_k);
    for j = 1:nb_data_in_prediction_media_k
      cell_sheet9(1,kl_index) = {"Prediction"};
      cell_sheet9(2,kl_index) = media(index_prediction_media(k));
      cell_sheet9(3,kl_index) = all_days_for_media_k(j);
      kl_index =  kl_index + 1;
    end
  end
  writecell({"Macroscopic rates from kinetic model"},name_excel_saving_file,'Sheet',9,'Range','A1');
  writecell(cell_sheet9,name_excel_saving_file,'Sheet',9,'Range','B1');
  writecell(num2cell((1:n_macro_rates)'),name_excel_saving_file,'Sheet',9,'Range','A4');
  writecell(num2cell([w_model_train,w_model_predict]),name_excel_saving_file,'Sheet',9,'Range','B4');

  % we save the cell specific rate computed from the colmun generation 
  cell_sheet10 = cell(3+n_macro_rates,N);
  kl_index = 1;
  for k = 1:length(nb_data_per_training_media)
    nb_data_in_training_media_k = nb_data_per_training_media(k);  
    kk = index_training_media(k);
    ind_days_media_k = (ic == kk);
    all_days_for_media_k = days_column(ind_days_media_k);
    for j = 1:nb_data_in_training_media_k
      cell_sheet10(1,kl_index) = {"Training"};
      cell_sheet10(2,kl_index) = media(index_training_media(k));
      cell_sheet10(3,kl_index) = all_days_for_media_k(j);
      kl_index =  kl_index + 1;
    end
  end
  for k = 1:length(nb_data_per_prediction_media)
    nb_data_in_prediction_media_k = nb_data_per_prediction_media(k);  
    kk = index_prediction_media(k);
    ind_days_media_k = (ic == kk);
    all_days_for_media_k = days_column(ind_days_media_k);
    for j = 1:nb_data_in_prediction_media_k
      cell_sheet10(1,kl_index) = {"Prediction"};
      cell_sheet10(2,kl_index) = media(index_prediction_media(k));
      cell_sheet10(3,kl_index) = all_days_for_media_k(j);
      kl_index =  kl_index + 1;
    end
  end
  writecell({"Macroscopic rates rates from CG"},name_excel_saving_file,'Sheet',10,'Range','A1');
  writecell(cell_sheet10,name_excel_saving_file,'Sheet',10,'Range','B1');
  writecell(num2cell((1:n_macro_rates)'),name_excel_saving_file,'Sheet',10,'Range','A4');
  writecell(num2cell([Macro_rates_train,Macro_rates_predict]),name_excel_saving_file,'Sheet',10,'Range','B4');
