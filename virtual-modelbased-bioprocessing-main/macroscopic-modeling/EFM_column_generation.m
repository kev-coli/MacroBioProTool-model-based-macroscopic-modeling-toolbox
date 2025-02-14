%   [Eb,Wb,cost_MP] = EFM_column_generation(Y,K,Eb,Aint,Dirrev,tolerance_EFM,use_mosek) 
%
%   This function is the column generation algorithm applied to macroscopic modeling
%   It computes the subset of EFMs Eb (one column = one EFM) such it minimizes
%   the cost norm(Y-K*Eb*Wb,2)^2 where Wb is the matrix which contains the
%   macro-rates, Y is the data matrix and K is the matrix involved in teh
%   column generation problem (for instance it can be Ameas, the measures stoichiometric matrix)
%
%   @inputs:
%       - Y : data to be fitted with the column generation
%       - K : extracellular stoichiometric matrix of the measured metabolites
%       - Eb : matrix containing the EFM which have already been computed
%       - Aint : intracellular stoichiometric matrix
%       - Dirrev : matrix used to impose irreversibility constraints on the fluxes
%       - use_mosek : boolean. If equal to 1, the column generation algorithm will use MOSEK. IT MUST BE INSTALLED AND THE PATHWAY INDICATED IN MAIN_MODELING.
%       - tolerance_EFM : threshold under which the modeling error for the EFM computation is judged negligible enough : all the most relevant EFMs are therefore computed
%       - text_to_print : text printed in the command window during the column generation algorithm
% 
%   @outputs:
%       - Eb : EFM matrix obtained after the column generation
%       - Wb : matrix containing the final computed macroscopic rates
%       - cost_MP : squared error cost of the EFM model with the data in Y

function [Eb,Wb,cost_MP] = EFM_column_generation(Y,K,Eb,Aint,Dirrev,tolerance_EFM,use_mosek) 

  % Column generation
  optimum_reached_with_column_generation = 0; % variable which tells the algorihm if the minimum has been reached
  while(optimum_reached_with_column_generation == 0) 
    clc
    fprintf("Number of generated EFMs : %d.\n",size(Eb,2))  
    [Enew,cost_SP,cost_MP,Wb] = master_and_subproblem(Y,K,Eb,Aint,Dirrev,use_mosek); % Solve the master and subproblems
    if(cost_SP < 0 && cost_MP > tolerance_EFM)  % If satisfied = minimum not reached
      Eb = [Eb,Enew]; % Concatenate the newly compute EFMto the old set EB of EFMs
    else
      optimum_reached_with_column_generation = 1; % The minimum has been reached, this variable becomes equal to 1 to stop the while loop.
    end
  end

end