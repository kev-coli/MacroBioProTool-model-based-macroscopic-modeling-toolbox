

use_mosek = 0; % if use_mosek is set to 1, it will use MOSEK, otherwise it should be set to 0.                 
               % BE CAREFUL : MOSEK 10.0 of after needs to be installed in
               % order to use it. Check their website : https://docs.mosek.com/latest/toolbox/install-interface.html
               % THE PATHWAY TO MOSEK FOLDER WHERE mosekdiag.m IS PRESENT ALSO NEEDS TO BE SPECIFIED AT THE
               % BEGINNING OF THIS .m FILE (see line 3)

%%%%%%%%% IMPORTANT FOR THE USE OF MOSEK !!!!!!
%%%%%%%%% ADD HERE THE PATHWAY TO MOSEK TOOLBOX
addpath('C:\Program Files\Mosek\10.2\toolbox\r2017aom')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Description of the code %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code allows to compute a macroscopic model of the
% cell kinetics using the data of concentration c_ext of 
% some extracellular metabolites. The model is such that 
%       q_ext = Aext*E*w(c_ext,theta)
% with
%   - q_ext: the cell specific rate data of the measured
%            metabolites (must be in the file)
%   - A_ext: the extracellular stoichiometric matrix
%   - E: the set of elementary flux modes (EFMs) computed during Step 2
%   - w: the vector of macroscopic rates identified as Monod functions
%   - theta: the parameter vector containing all the activation, inhibition and maximal rate constants to be identified     

% The identification method can be divided into 4 steps :
% STEP 1 : the data in the excel files or mat files (metabolic network, cext and qext data) are treated 
% STEP 2 : a subset of EFMs E is computed using the data qext and the metabolic network.
%          Additionnally, the macroscropic rate data are also computed.   
% STEP 3 : theta is identified from the macroscopic rate data obtained during Step 2 and the cext data
%          For this step, the code first involves a model selection (Step
%          3.a) by using l1 norm regularization.
%          parameter tuning with Bayesian estimation (Step 3.b) 
%          and modulation function reduction (Step 3.c)
% STEP 4 : treatment of the results such as plots, error computation and
%          saving important files for model-based optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% FILES REQUIREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each of these steps require three files whose pathway need to be specified by the user

% 1) the file containing the qext data (.xlsx or .mat file)
%       .xlsx file : Except for the two first columns, each column contains the cell specific rqtes qext for a given metabolite
%                    Cells C1-XXX1 must contain the metabolites name (SAME name as in the stoichiometric matrix),
%                    The first column (A) has the name of each experimental condition
%                    The second column (B) contains the measurement time (in days) for each experimental condition
%                    The rest (starting from C3) contains the qext values
%                    See fake data example for illustration
%       . mat file : it must contain three variables
%                       - qext of array type : contains the values of the cell specific rates. Each column represents a metabolite, the rows, the different measurements
%                       - 

% 2) the file containing the cext data (.xlsx or .mat file)
%    SAME STRUCTURE FOR BOTH THE EXCEL AND MATLAB FILES AS THE FILE CONTAINING THE qext DATA

% 3) the metabolic network file containing the stoichiometric matrix, reactions reversibility and the intracellular or extracellular property of each metabolite


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% VARIABLES FOR STEP 1 - DATA TREAMTENT %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VARIABLES FOR STEP 1 - FILE DIRECTORY
% directory_X variables must be in char or string!
% They must be either .xlsx ot .mat files!
% WARNING: see the example files in order to have the right format

% Directory for the file contaning the qext data
directory_file_cell_specific_rate_data_qext = '../data/K-Net data/fake_data_example_qext2.xlsx';

% Directory for the file contaning the cext data
directory_file_concentration_data_cext = '../data/K-Net data/fake_data_example_cext.xlsx';

% Directory for the file contaning the metabolic network
directory_file_metabolic_network = '../data/K-Net data/fake_data_example_metabolic_network.xlsx';

prediction_media = {'C6'}; %  Put here in char the conditions which are excluded from the model training stage, they will only be used for prediction

% Smoothing : it is possible to smooth the  data  cext and qext per conditon
% The code uses the Matlab function smooth
smoothing = 0; % Smoothing of the data cext and qext -> 1 : yes, 0 : no
coeff_smoothing = 0.05; % Span for the smoothing (see smooth function on matlab for more details)

% Normalization of qext 
normalization = 0; % Normalization of the qext data -> 1 : yes, 0 : no
normalization_matrix = []; % Weighting square matrix used for the normalization of the qext data.
                           % It must have the same dimension as the number
                           % of extracellular metabolites present in the
                           % qext data file
                           % If normalization = 1 and normalization_matrix = [], the code will automatically normalize
                           % with the average value of the data

average_data_per_condition = 0; % This computes the average of all the data of a given condition. Interesting for data at steady-state. Set it to 1 for averaging the data, 0 otherwise.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% VARIABLES FOR STEP 2 - EFM COMPUTATION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The code can use MOSEK optimization tool which is faster and more efficient 

computation_method_EFM = 'sequential'; % 3 possible choices : 'global', 'union' or 'sequential'
tolerance_EFM = 1e-10; % Tolerance for the computation of the EFMs
reduction_of_EFM = 1; % 2 possible choices : 1 for reduction, 0 without reduction
factor_error_reduction = 10; % It corresponds to the PERCENTAGE of cost degradation allowed for the EFM reduction.
                             % If cost_MP is the least-squares error after
                             % column generation, the code will remove EFMs
                             % until the least-squares error is above a
                             % threshold cost given by
                             % (100 + fact_error_reduction)/100*cost_MP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% VARIABLES FOR STEP 3 - MONOD KINETICS IDENTIFICATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mets_not_included_in_Monod_kinetics = {'X'}; % Put here in char the name of the metabolites which are not included in the Monod kinetics expression

% As aformentioned, Step 3 follows three substeps (model selection in Step 3.a, Bayesian estimation in Step 3.b and model reduction in Step 3.c)
% Several variables need to be defined for each of these steps.

%% Variables for Step 3.a: model selection with regularization %%

% In step 3.a, the code solves 
%    min ||w(k) - w(cext(k),theta)||^2 + lambda*(sum of activation and inhibition parameters)
% lambda is the regularization parameter (non-negative)
% Regularization allows to set a lot of kinetic parameters to 0.

regularization = 1;  % To add regularization, set regularization to 1. Otherwise set it to 0. 
                     % WARNING : if there are more parameters in the Monod
                     % model than data, regularization is essential !
lambda_grid = [0,0.01]; % Specify the grid for the regularization parameter 
number_multistart_points = 3; % number of multistart points for the regularization optimization problem
                              % the higher, the better for global minimum
                              % computation. However the code will take
                              % more time.

%% Variables for Step 3.b: model selection with regularization %%

% In step 3.b, Bayesian estimation is used in order to re-tune the parameters
% which have not been set to 0 after the model selection step (Step 3.a).
% This can possibly improve the model accuracy since regularization can
% introduce a large bias 

number_max_of_iterations_EM = 100;  % Number of maximal iterations for the Expectation Maximization algorithm
burn_in_period_EM = 10000; % Burn-in period for the Metropolis-Hastings algorithm (it should be large enough for a proper warm-up of the algorithm)
number_samples_per_iteration_EM = 1000; % Number of samples for the hyperparemeter updates  (the higher the better for parameter accuracy)
number_trials_sampling_EM = 1; % Number of trials for the Metropolis Hastings algorithm  (the higher the better for parameter convergence)
perturbation_variance_EM = 1; % perturbation of the variance hyperparameters for more exploration of the parameter space (advice : do not exceed 0.5 otherwise the convergence might not happen)

%% Variables for Step 3.c: modulation function reduction %%

% In step 3.c, each identified modulation function is analyzed in order to
% check if they can be reduced in complexity. Two tests are performed for
% each modulation function

% Test 1 : we first check if it can be a neutral effect by checking if
% the minimal and maximal values of the modulation function computed with
% the cext data have a large relative error. If it is below a given
% threshold defined by the user (with the variable
% variation_tolerance_neutral_modulation_function), then the modulation function is 
% considered as neutral effect and both the identified activation and inhibition 
% parameters are set to 0.

% Test 2 : if Test 1 was not successful, we check if the modulation
% function can be reduced into an inhibition or activation effect. From the
% output of the identified modulation function after Step 3.b with the cext data, we
% compute the best activation and inhibition model (computed with the same cext data)
% which minimizes the squared error. We take the model with the largest fit
% (in percentage). If this largest fit is above a percentage threshold
% chosen by the user (threshold_fit_activation_or_inhibition), then the modulation function is 
% reduced into the corresponding effect. 

threshold_variation_neutral_effect = 5; % Variation threshold percentage for the reduction of the double-compnent kinetics into neutral effect (Test 1)
threshold_fit_activation_or_inhibition = 95; % Threshold percentage for the reduction of the double-compnent kinetics into either activation of inhibition (Test 2 if Test 1 not successful)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% VARIABLES FOR STEP 4 - PLOT, SAVE and TREAT THE RESULTS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The user specifies here the name of the folder and the file (with the same name stored in the variable name_save_file) containing the save (in a .mat and a .xlsx file)
% IMPORTANT INFORMATION : the folder and the files within will be saved in the folder models of the toolbox

name_save_file = "K-Net";

% For the plots of the cell specific rates, the user can choose the number
% of plots to put per rows and columns
number_of_plots_per_columns = 2; % number of plots per columns
number_of_plots_per_rows = 2;  % number of plots per rows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           NO MORE VARIABLES TO BE DEFINED            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code runs STEPS 1, 2, 3 and 4
run_all_the_steps
