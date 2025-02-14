%   [data_qext, data_cext, data_stoich] = Step_1_file_treatment(all_file_directories,options_data)     
%
%   This code performs several functons :
%   - checks if the files (.mat, .xlsx or .xls) containing the cell specific rate
%     data qext, the concentration data cext and the metabolic network exist. The directories of these files must be defined
%     in MAIN_MODELING with the following corresponding variables: directory_file_cell_specific_rate_data_qext,
%     directory_file_concentration_data_cext and directory_file_metabolic_network
%   - loads/extracts necessary informations for the modeling (such as stroichiometric matrix, the number of data, reversibility of the reactions ,etc)
%   - smoothes the data cext and qext (or not) according to the choice of the user made in MAIN_MODELING.m
%   - normalizes the qext data (or not) according to the choice of the user made in MAIN_MODELING.m
%   - stores all the necessary information for modeling in the structures data_qext, data_cext and data_stoich
%
%   @inputs:
%       - all_file_directories : a cell vector of dimension 3 containing the directory of the files with the qext data, cext data and the metabolic network, in that order.
%       - options_data : a structure with three entries
%               -> options_data.smoothing : boolean. 0 if the user does not want to smooth the data cext and qext, 1 if he does
%               -> options_data.coeff_smoothing : span for the smoothing as used in the Matlab fonction smooth.
%               -> options_data.normalization : boolean. 0 if the user does not want to normalize the data qext, 1 if he does
%               -> options_data.normalization_matrix : normalization matrix defined by the user. If it is empty, the code will consider a diagonal matrix with the average of the data in the diagonal elements
%               -> options_data.average_data_per_condition : boolean. If it is equal to 1 (resp. 0), all the data of the same ocndition are averaged (resp. not averaged). 
%               -> options_data.prediction_media : cell vector containing the name of the conditions used for prediction
%
%   @outputs:
%       - data_qext : structure containing all the information about the data qext which is used in STEPs 2, 3 and 4
%               -> data_qext.ia : vector used in order to identify the condition to which each measurement corresponds to (combined with data_qext.ic) 
%               -> data_qext.ic : vector used in order to identify the condition to which each measurement corresponds to (combined with data_qext.ia)  
%               -> data_qext.media : cell vector with the name of each experimental condition identified in the qext qnd cext data files
%               -> data_qext.n_media : number of different experimental conditions
%               -> data_qext.index_training_media : vector gathering all the indices of the conditions in media which are dedicated for model training
%               -> data_qext.index_prediction_media : vector gathering all the indices of the conditions in media which are dedicated for model prediction
%               -> data_qext.N : number of measurements (data)
%               -> data_qext.mets_meas_in_qext : cell with all the names of the extracellular metabolites whose cell specific rate is measured
%               -> data_qext.qext : matrix of the cell specific rate data. Each row corresponds to a given metabolite and the columns correspond to the measurement time instants
%               -> data_qext.qext_train : matrix of the cell specific rate TRAINING data (averaged, smoothed and/or normalized)
%               -> data_qext.qext_predict : matrix of the cell specific rate PREDICTION data (averaged, smoothed and/or normalized)
%               -> data_qext.media_column : cell vector containing the name of the media/condition to which each measurement comes from
%               -> data_qext.days_column : cell vector containing the measurement time instants of each condition
%       - data_cext : structure containing all the information about the data cext which is used in STEPs 2, 3 and 4
%               -> data_cext.ic : same as data_qext.ic
%               -> data_cext.ic_training : vector used in order to identify the TRAINING condition to which each measurement corresponds to 
%               -> data_qext.mets_meas_in_cext : cell with all the names of the extracellular metabolites whose concentration is measured
%               -> data_qext.media : cell vector with the name of each experimental condition identified in the qext qnd cext data files
%               -> data_qext.cext : cell vector with the name of each experimental condition identified in the qext qnd cext data files
%               -> data_cext.cext_train : matrix of the concentration TRAINING data (averaged, smoothed and/or normalized)
%               -> data_cext.cext_predict : matrix of the concentration PREDICTION data (averaged, smoothed and/or normalized)
%       - data_stoich : structure containing all the information about the metabolic network which is used in STEPs 2, 3 and 4
%               -> data_stoich.mets : cell vector with all the metabolites involved in the metabolic network   
%               -> data_stoich.A : stoichiometric matrix. The rows correspond to the metabolites and the columns to the reactions 
%               -> data_stoich.n_v : number of reactions. It also corresponds to the number of columns in the stoichiometric matrix A
%               -> data_stoich.is_ext : boolean vector of same dimension as mets. It the j-th entry of is_ext is 1 (resp. 0), then the j-th metabolite in mets is extracellular (resp. intracellular).
%               -> data_stoich.mets_ext : cell vector which gathers the extracellular metabolites 
%               -> data_stoich.mets_int : cell vector which gathers the intracellular metabolites 
%               -> data_stoich.n_mets_int : number of intracellumar metabolites
%               -> data_stoich.rev_bool : boolean vector whose dimension is equal to the number of reactions (n_v). It the j-th entry of rev_bool is 1 (resp. 0), then the j-th reaction is reversible (resp. irreversible).
%               -> data_stoich.ind_rev : vector made up of the indices of the reversible reactions in the flux vector 
%               -> data_stoich.ind_irrev : vector made up of the indices of the irreversible reactions in the flux vector 
%               -> data_stoich.n_v_rev : number of reversible reactions (equal to the dimension of ind_rev)
%               -> data_stoich.n_v_irrev : number of irreversible reactions (equal to the dimension of ind_irrev)
%               -> data_stoich.Dirrev : matrix of dimension n_v_irrev*n_v such that v_irrev = D_irrev*v with v_irrev the subvector of the flux vector v containing only the irreversible fluxes. Useful in order to implement the constraint Dirrev*v >=0. 
%               -> data_stoich.Aext : extracellular stoichiometric matrix. The rows correspond to the EXTRACELLULAR metabolites and the columns to the reactions
%               -> data_stoich.Aint : intracellular stoichiometric matrix. The rows correspond to the  INTRACELLULAR metabolites and the columns to the reactions
%               -> data_stoich.is_meas : boolean vector whose dimension equals the number of extracellular metabolites. If the j-th entry of is_meas is 1 (resp. 0), then the cell specific rate of the j-th extracellular metabolite in mets_ext is measured (resp. not measured)
%               -> data_stoich.Ameas : extracellular stoichiometric matrix with only the measured extrecallular metabolites.
%               -> data_stoich.n_meas : number of extrecallular metabolites whose cell specific rate is measured
%               -> data_stoich.mets_meas_in_stoich : measured metabolites which could be identified in the metabolic network

function [data_qext,data_cext,data_stoich] = Step_1_file_treatment(all_directories,options_data)
  
  clc
  fprintf("Step 1: Treatment of all the files.\n")
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Treatment of the file containing the extracellular cell specific rate data qext %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  directory_file_qext = char(all_directories{1});
  if(exist(directory_file_qext) == 0) % Test to check if the file exists in the directory specified by the user.
    error("The file for the extracellular cell specific rate data qext does not exist.\n") % If not, the file does not exist -> the code stops here
  end
  [~,~,extension_file_qext] = fileparts(directory_file_qext); % If yes, get the extension of the file
  
  fprintf("File containing the data qext...") 
  pause(1)

  % Loading the file containing qext data
  % It will also store necessary variables in the structure data_qext
  if(isequal(extension_file_qext,'.xlsx') || isequal(extension_file_qext,'.xls')) % If the file is from excel
    content_file_qext = readtable(directory_file_qext); % read the excel file
    headers = content_file_qext.Properties.VariableNames; % Get the first row from the excel
    length_headers = length(headers); % Compute the number of columns
    mets_meas_in_qext = headers([3:length_headers]); % Get the measured metabolic names from xlsx files
    media_column = content_file_qext(1:end,1); % Get the name of the experimental conditions
    days_column = content_file_qext(1:end,2); % Get the name of the experimental conditions
    [media, ia, ic] = unique(media_column.(1),'stable'); % Extract the number of media + their names
                                                         % ia: first index where a given condition appears first
                                                         % ic: indices where a given condition appear in the data  
    media_column_qext = media_column.(1); 
    days_column_qext = days_column.(1); 
    data_qext.media_column = media_column_qext;
    data_qext.days_column = days_column_qext;
    qext = table2array(content_file_qext(:,[3:length_headers])); % uptake/secretion rates
    data_qext.ia = ia; % store ia in data_qext
    data_qext.ic = ic; % store ic in data_qext
    data_cext.ic = ic; % store ic in data_cext
    data_qext.media = media;
    n_media = length(media); 
    data_qext.n_media = n_media;
    data_qext.mets_meas_in_qext = mets_meas_in_qext; % store mets_meas_in_qext in data_qext  
    
    % Treatment analysis between training and prediction media
    prediction_media = options_data.prediction_media; % load the names of the prediction media specified by the user
    [~,ind_prediction_media] = ismember(prediction_media,media); % identify the names in the cell vector media which are the same as in prediction_media
    index_prediction_media = nonzeros(ind_prediction_media)'; % indices of the media in the cell vector media which are dedicated for prediction 
    index_training_media = setdiff(1:n_media,ind_prediction_media); % indices of the media in the cell vector media which are dedicated for training
    data_qext.index_training_media = index_training_media;
    data_qext.index_prediction_media = index_prediction_media;
  else
    error("The file does not have the right format. The format should be .mat, .xls or .xlsx.")
  end
  
  % Treatment of possible missing data                                                                                     
  [rows_q,cols_q,~] = find(isnan(qext)); % Look for the (i,j)-entries of qext with entries equal to NaN (missing data)
  %  Fill missing values in qext by replacing by the average of the data
  if ~isempty(rows_q) % If we have some missing data in qext (some NaN entries)
    for kk = 1:length(rows_q) % For each missing data (NaN entry)
      tempc = qext(ic == ic(rows_q(kk)),cols_q(kk)); % We take all the samples of the medium which has the kk-th missing data
      qext(rows_q(kk),cols_q(kk))= mean(tempc(~isnan(tempc))); % We take the mean of the samples of this medium (excluding NaN entries) and we replace the missing data by this value
    end
  end
  N = size(qext,1); % Total number of measurements
  data_qext.N = N; % store N in data_qext
  fprintf("OK!\n")
  pause(1)

  %% Treatment of the file containing the extracellular concentration data cext
  directory_file_cext = char(all_directories{2});
  if(~exist(directory_file_cext)) % Test to check if the file exists in the directory specified by the user.
    error("The file for the extracellular concentrations does not exist.\n") % If not, the file does not exist -> the code stops here
  end
  [~,~,extension_file_cext] = fileparts(directory_file_cext); % If yes, get the file  extension   

  fprintf("File containing the data cext...") 
  pause(1)

  if(isequal(extension_file_cext,'.xlsx') || isequal(extension_file_cext,'.xls')) % If the file is from excel, load its content
    content_file_cext = readtable(directory_file_cext);
    headers = content_file_cext.Properties.VariableNames; % Get the first row from the excel file
    length_headers = length(headers); % Compute the number of columns
    media_column = content_file_cext(1:end,1); % Get the name of the experimental conditions
    days_column = content_file_cext(1:end,2); % Get the name of the experimental conditions
    media_column_cext = media_column.(1); 
    days_column_cext = days_column.(1); 
    if(~isequal(media_column_cext,media_column_qext) || ~isequal(days_column_cext,days_column_qext)) % We check here if both first columns of the cext and qext files are the same (i.e., same measurement time instants and condition names)
      error("The files contaning the cext and qext data must have the same measurement time instants and condition names.\n This requires that the both first columns of both excel files are the same between both files.")
    end
    mets_meas_in_cext = headers([3:length_headers]); % Get the name of the extracellular metabolites whose concentrations are measured
    cext = table2array(content_file_cext(:,[3:length_headers])); % uptake/secretion rates
  else
    error("The file does not have the right format. The format should be .xls or .xlsx.")
  end
  
  % Treatment of possible missing data
  [rows_c,cols_c,~] = find(isnan(cext)); % Look for the (i,j)-entries of cext with entries equal to NaN (missing data)
  %  Fill missing values in cext by replacing by the average of the data
  if ~isempty(rows_c) % If we have some missing data in cext (some NaN entries)
    for kk = 1:length(rows_c) % For each missing data (NaN entry)
      tempc = cext(ic == ic(rows_c(kk)),cols_c(kk)); % We take all the samples of the medium which has the kk-th missing data
      cext(rows_c(kk),cols_c(kk))= mean(tempc(~isnan(tempc))); % We take the mean of the samples of this medium (excluding NaN entries) and we replace the missing data by this value
    end
  end
  fprintf("OK!\n")
  pause(1)

  %% File of the metabolic network
  % Does the file exist? 
  directory_file_metabolic_network = char(all_directories{3});
  if(~exist(directory_file_metabolic_network))
    error("The file for the metabolic network does not exist.\n")
  end
  [~,~,extension_file_metabolic_network] = fileparts(directory_file_metabolic_network); 
  
  fprintf("File containing the metabolic network...") 
  pause(1)

  % This code loads the corresponding file if it has a .mat, .xls or .xlsx format. 
  % In both cases, the content is stored in the variable content_file_metabolic_network (structure for .mat file, cell for .xls and .xlsx format)
  % The content of content_file_metabolic_network is then treated in order to extract:
  %   - the irreversibility of each reaction (boolean vector is_irrev). The j-th entry of this vector is 1 if the j-th reaction is irreversible, 0 if reversible.
  %   - the name of the metabolites (stored in the cell vector mets)
  %   - the nature of these metabolites (extracellular of intracellular), this informaion is stored in the boolean vector (is_ext). The j-th entry of this vector is 1 if the j-th metabolite is extracellular, 0 if intracellular.
  %   - the stoichiometric matrix A (double array) whose rows correspond to the metabolites and columns to the reactions 
  if(isequal(extension_file_metabolic_network,'.mat'))  % If the file is .mat, load the content
    content_file_metabolic_network = load(directory_file_metabolic_network);
    mets = content_file_metabolic_network.mets; % cell containing the name of the metabolites
    data_stoich.mets = mets; % storing mets in the structure data_stoich which will be used later
    is_ext = content_file_metabolic_network.is_ext; % boolean vector containing the extracellular and intracellular 
    data_stoich.is_ext = is_ext;
    A = content_file_metabolic_network.A;
    data_stoich.A = A;
    is_irrev = content_file_metabolic_network.is_irrev;
    data_stoich.is_irrev = is_irrev;
  elseif(isequal(extension_file_metabolic_network,'.xlsx') || isequal(extension_file_metabolic_network,'.xls')) % If the file is from excel, load the content
    content_file_metabolic_network = readtable(directory_file_metabolic_network);
    mets = table2array(content_file_metabolic_network(2:end,2)); % metabolites acronyms 
    data_stoich.mets = mets;
    extracellular_or_intracellular = table2array(content_file_metabolic_network(2:end,5)); % are the metabolites extracellular or intracellular ('ext' = extracellular, 'int' = intracellular)
    is_ext = (string(extracellular_or_intracellular) == "ext"); % boolean vector for extrecallular or intracellular metabolites  (1 = extracellular, 0 = intracellular)
    data_stoich.is_ext = is_ext;
    A = table2array(content_file_metabolic_network(2:end,12:end)); % stoichiometric matrix
    data_stoich.A = A;
    is_irrev = (table2array(content_file_metabolic_network(1,12:end)) == 1)';
  else
    error("The file does not have the right format. The format should be .mat, .xls or .xlsx.")
  end
  fprintf("OK!\n")
  pause(1)
  
  fprintf("Preparation of all the necessary variables for modeling...")
  pause(2)

  % Load the smoothing and normalization options chosen by the user
  smoothing = options_data.smoothing;
  coeff_smoothing = options_data.coeff_smoothing;
  normalization = options_data.normalization;
  normalization_matrix = options_data.normalization_matrix;
  average_data_per_condition = options_data.average_data_per_condition;

  % Split the metabolites between the extracellular ones and the intracellular ones
  % is_ext is a boolean vector. If its j-th entry is equal to 1, then the j-th metabolite 
  % in the vector mets is extracellular. If 0, it is intracellular.
  is_ext = data_stoich.is_ext;
  mets_ext = mets(is_ext); % Gather all extracellular metabolites in mets_ext thanks to the index vector is_ext
  data_stoich.mets_ext = mets_ext;
  mets_int = mets(~is_ext); % Gather all extracellular metabolites in mets_int
  data_stoich.mets_int = mets_int;

  % Look for the reversible reactions and count them
  % is_irrev is a boolean vector. If its j-th entry is equal to 1, then the j-th reaction 
  % is irreversible. If 0, it is reversible. 
  n_v = length(is_irrev); %  Number of fluxes (= number of columns in the stoichiometric matrix)
  data_stoich.n_v = n_v;
  n_v_irrev = sum(is_irrev == 1); %  Number of irreversible fluxes 
  data_stoich.n_v_irrev = n_v_irrev;
  rev_bool = ~is_irrev; % rev_bool is a boolean vector. If its j-th entry is equal to 1, then the j-th reaction 
                          % is reversible. If 1, it is irreversible. 
  data_stoich.rev_bool = rev_bool;                        
  n_v_rev = sum(rev_bool == 1); %  Number of reversible fluxes  
  data_stoich.n_v_rev = n_v_rev;  
  ind_irrev = nonzeros(diag(1:n_v)*is_irrev); % Indices of the ractions which are irreversible 
  data_stoich.ind_irrev = ind_irrev; 
  ind_rev = nonzeros(diag(1:n_v)*(~is_irrev)); % Indices of the ractions which are reversible  
  data_stoich.ind_rev = ind_rev; 
  % Construction of the matrix Dirrev which is such that v_irrev = Dirrev*v
  % where v is the vector contaning the fluxes of all the reactions and v_irrev is the subvector of v containing
  % only the fluxes of the irreversible reactions
  Dirrev = diag(is_irrev); Dirrev(~any(Dirrev,2),:) = [];
  data_stoich.Dirrev = Dirrev; 

  data_stoich.A = A; 
  Aext = A(is_ext,:); % Extracellular stoichiometric matrix 
  data_stoich.Aext = Aext; 
  Aint = A(~is_ext,:); % Intracellular stoichiometric matrix
  data_stoich.Aint = Aint; 
  n_mets_int = size(Aint,1); % Number of intracellular metabolites
  data_stoich.n_mets_int = n_mets_int; 

  % We look for the extracellular metabolites which are measured and we
  % look for their indices in the stoichiometric matrix Aext
  [is_meas,~] = ismember(mets_ext,mets_meas_in_qext); % We look for the extracellular metabolites which are measured
  data_stoich.is_meas = is_meas; 
  mets_meas_in_stoich = mets_ext(is_meas); % Extract the MEASURED extracellular metabolites from teh stoichiometric matrix
  data_stoich.mets_meas_in_stoich = mets_meas_in_stoich; 
  Ameas = Aext(is_meas,:); % Stoichiometric matrix of the MEASURED extracellular metabolites
  data_stoich.Ameas = Ameas; 
  n_meas = size(Ameas,1); % Number of MEASURED extracellular metabolites
  data_stoich.n_meas = n_meas; 

  % Rearrange the measured metabolites in the same order of appearance in the measured stoichiometric matrix Ameas
  % To do so, the function ismember(A,B) can also return an array locb containing the
  % lowest absolute index in B for each element in A
  % Hence, the order of the metabolites in A is kept
  [~, locb] = ismember(mets_meas_in_stoich,mets_meas_in_qext); 
  qext = qext(:,locb); % Right order of qext for the MEASURED metabolites
  data_qext.qext = qext;
  mets_meas_in_qext = mets_meas_in_qext(locb);
  data_qext.mets_meas_in_qext = mets_meas_in_qext;

  % Rearrange the measured metabolites in cext in the same order of
  % appearance in the measured stoichiometric matrix Ameas (same order as the rows of Ameas)
  [~, loca] = ismember(mets_meas_in_stoich,mets_meas_in_cext);
  loca = nonzeros(loca);
  cext = cext(:,loca); % Right order of the columns of cext (same as for the rows of the stoichiometric matrix Ameas)
  mets_meas_in_cext = mets_meas_in_cext(loca); 
  data_cext.cext = cext;  
  data_cext.mets_meas_in_cext = mets_meas_in_cext;
  data_cext.media = media;
  
  % Identification the number of measurement time instants for each training and prediction medium
  nb_data_per_training_media = []; % vector of same size as the number of training conditions which indicate the number of measurements for each training condition
  nb_data_per_prediction_media = []; % vector of same size as the number of prediction conditions which indicate the number of measurements for each prediction condition
  ic_training = []; kl_train = 1;
  if(~isempty(index_prediction_media))
    for k = index_prediction_media
      nb_data_per_prediction_media = [nb_data_per_prediction_media,sum(ic == k)];
    end
  end
  for k = index_training_media
    ind_k = (ic == k);
    nb_data_per_training_media = [nb_data_per_training_media,sum(ind_k)];
    ic_training = [ic_training;kl_train*ones(sum(ind_k),1)];
    kl_train = kl_train + 1;
  end
  data_cext.ic_training = ic_training;
  data_qext.nb_data_per_training_media = nb_data_per_training_media;
  data_qext.nb_data_per_prediction_media = nb_data_per_prediction_media;

%%%%%% Smoothing, normalization and averaging
  qext_treated = qext; % qext_treated is equal to qext with smoothing, averaged and normalization (depending of the choice of the user)
  cext_treated = cext; % qext_treated is equal to qext with smoothing, averaged and normalization (depending of the choice of the user)
    
  % Normalization of the data qext
  if(normalization)
    if(isempty(normalization_matrix))    
      normalization_matrix = abs(diag(mean(qext,1))^-1);  
    end
    qext_treated = qext_treated*normalization_matrix; 
  end

  % Smoothing or averaging (per condition) of the data qext
  if(smoothing == 1 && average_data_per_condition == 0)
    for ii = 1:size(qext_treated,2)
      qext_treated(:,ii) = smooth(qext(:,ii),coeff_smoothing,'loess');
    end
  elseif(smoothing == 0 && average_data_per_condition == 1)
    qext_averaged = [];
    for ii = 1:n_media
      qext_averaged = [qext_averaged;mean(qext(ic == ii,:),1)];
    end
    qext_treated = qext_averaged;
    qext_train = qext_treated(index_training_media,:);
    if(~isempty(index_prediction_media))
      qext_predict = qext_treated(index_prediction_media,:);
    else
      qext_predict = [];  
    end
  elseif(smoothing == 0 && average_data_per_condition == 0)
    qext_train = []; qext_predict = [];
    if(~isempty(index_prediction_media))
      for k = index_prediction_media
        ind_k = (ic == k);
        qext_predict = [qext_predict;qext_treated(ind_k,:)];
      end
    end
    for k = index_training_media
      ind_k = (ic == k); 
      qext_train = [qext_train;qext_treated(ind_k,:)];
    end    
  end  
  data_qext.qext_train = qext_train;
  data_qext.qext_predict = qext_predict;
 
  % Smoothing or averaging of the data cext
  % + split between training and prediction data
  if(smoothing == 1 && average_data_per_condition == 0)
    for ii = 1:size(cext,2)
      cext_treated(:,ii) = smooth(cext(:,ii),coeff_smoothing,'loess');
    end
  elseif(smoothing == 0 && average_data_per_condition == 1)
    cext_averaged = [];
    for ii = 1:n_media
      cext_averaged = [cext_averaged;mean(cext(ic == ii,:),1)];
    end
    cext_treated = cext_averaged;
    cext_train = cext_treated(index_training_media,:);
    if(~isempty(index_prediction_media))
      cext_predict = cext_treated(index_prediction_media,:);
    else
      cext_predict = [];  
    end
  elseif(smoothing == 0 && average_data_per_condition == 0)
    cext_train = []; cext_predict = [];
    if(~isempty(index_prediction_media))
      for k = index_prediction_media
        ind_k = (ic == k);
        cext_predict = [cext_predict;cext_treated(ind_k,:)];
      end
    end
    for k = index_training_media
      ind_k = (ic == k); 
      cext_train = [cext_train;cext_treated(ind_k,:)];
    end  
  end
  data_cext.cext_train = cext_train;
  data_cext.cext_predict = cext_predict;
  
  fprintf("Done!\n")
  pause(1)
  fprintf("Step 2 with EFM computation is about to start.\n")
  pause(3)

  
end