%   [Enew,cost_SP,cost_MP,Wb] = master_and_subproblem(Y,K,Eb,Aint,Dirrev,use_mosek)
%
%   This code solves the masterproblem (macro-rate computation) and the subproblem (computation of a new EFM candidate) for the coluln generation
%   For the masterproblem, MOSEK can be used. WARNING : IT MUST BE INSTALLED
%
%   @inputs:
%       - Y : data to be fitted with the column generation
%       - K : extracellular stoichiometric matrix of the measured metabolites
%       - Eb : matrix containing the EFMd which have already been computed
%       - Aint : intracellular stoichiometric matrix
%       - Dirrev : matrix used to impose irreversibility constraints on the fluxes
%       - use_mosek : boolean. If equal to 1, the column generation algorithm will use MOSEK. IT MUST BE INSTALLED AND THE PATHWAY INDICATED IN MAIN_MODELING.
% 
%   @outputs:   
%       - Enew : new EFM candidate
%       - cost_SP : subproblem cost
%       - cost_MP : masterproblem cost (data fit)
%       - Wb : matrix containing the macroscopic rates computed after the masterproblem

function [Enew,cost_SP,cost_MP,Wb] = master_and_subproblem(Y,K,Eb,Aint,Dirrev,use_mosek)
    
      options_qp = optimoptions('quadprog','Display','none','ConstraintTolerance',1e-10,'MaxIterations',10^6);
      options_lin = optimoptions('linprog','Display','none');
      options_lin.Preprocess = 'none';
      N = size(Y,2); % number of data in Y
      [n_mets_int,n_v] = size(Aint);
      n_Eb = size(Eb,2);  % Number of macroscopic reaction 
      n_v_irrev = size(Dirrev,1); % nulber of irreversible reactions
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%% Masterproblem %%%%%%%%%%%%%%%%%%%%%%%%%% 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % The master problem computes the macroscopic rates of the EFM which
      % have already been computed (in the matrix Eb)
      % It is a convex quadratic optimization to be solved with
      %      - the matlab function quadprog if use_mosek is set to 0 (not recommended)
      %      - the MOSEK second order cone programming optimizer if use_mosek is set to 1

      L = K*Eb; 
      H = L'*L + 1e-15*eye(n_Eb); % Quadratic matrix
      H = (H+H')/2;
      Wb = zeros(n_Eb,N); % matrix of macroscopic rate data
      cost_MP = 0; % masterproblem cost
      err = zeros(size(Y,1),N); % Error between the data and the EFM model Ameas*Eb*Wb
      parfor kk = 1:N
        Yk = Y(:,kk);  
        fk = -L'*Yk;  
        if(use_mosek)
          x = quadprog_mosek_as_second_order_programming(L,Yk,[],[],[],[],zeros(n_Eb,1),Inf*ones(n_Eb,1));  
        else
          x = quadprog(H,fk,[],[],[],[],zeros(n_Eb,1),[],[],options_qp);
        end
        x(x < 0) = 0;
        Wb(:,kk) = x;
        cost_MP = cost_MP + norm(Y(:,kk)-L*Wb(:,kk),2)^2;
        err(:,kk) = (L*Wb(:,kk) - Y(:,kk));
      end
      f0_all = K'*err;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%% Subproblem %%%%%%%%%%%%%%%%%%%%%%%%%%% 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % This part of the code solves the subproblem who has two purposes :
      %     - compute a new EFM candidate
      %     - check if the set of already computed EFM is enough for optimal data fit

      fo = sum(f0_all,2);
      fo = sparse(fo);
      A_lin = [-sparse(eye(n_v)),sparse(eye(n_v));
               sparse(eye(n_v)),sparse(eye(n_v));
               Dirrev,sparse(n_v_irrev,n_v)]; % Matrix implementing the linear inequality constraints
      b_lin = [sparse(2*n_v+n_v_irrev,1)];
      A_eq = [Aint,sparse(n_mets_int,n_v);sparse(1,n_v),ones(1,n_v)];  % Matrix implementing the linear equality constraints
      b_eq = [sparse(n_mets_int,1);1];
      f = [fo;sparse(n_v,1)]; % Vector of the linear objective function
      e_opt_new = linprog(f,-A_lin,-b_lin,A_eq,b_eq,[],[],options_lin); % compute the new EFM candidate
      [e_opt_new2] = improve_EFM_computation(e_opt_new,A_lin,b_lin,A_eq,b_eq,f);  % improve iteratively the accuracy of the EFM
      while(sum((e_opt_new2) > 0) < sum((e_opt_new) > 0))
        e_opt_new = e_opt_new2; 
        [e_opt_new2] = improve_EFM_computation(e_opt_new,A_lin,b_lin,A_eq,b_eq,f);
      end
      cost_SP = f'*e_opt_new2; % subproblem cost. If it is equal to 0, optimality has been reached already
      Enew = e_opt_new2(1:n_v,1); % new EFM candidate

      % This part of the code verifies if the new EFM candidate satisfied several properties to be considered an EFM
      if(abs(norm(Enew,1) - 1) < 1e-3) % First we verify that the computed EFM has a norm l1 equal to 1
        for i = 1:n_Eb % Then, we verify if the new EFM has not already been computed
          if(sum(abs(Enew - Eb(:,i)) <  1e-3) == n_v)  % If it is the case, we set cost_SP to 0 which will stop the coluln generation algorithm
            cost_SP = 0;
            Enew = zeros(n_v,1); 
            break         
          end
        end  
      else
        cost_SP = 0;
        Enew = zeros(n_v,1);  
      end     

end

function [e_opt_new] = improve_EFM_computation(e_opt,A_lin,b_lin,A_eq,b_eq,f)

% The purpose of this code is to compute again an EFM where all the entries which are equal to 0 are removed
% The idea is to iteratively reduce the EFM computation complexity (reduce the number of variables)
% It will increase the accuracy of the new optimum
%
%   @inputs:  
%       - e_opt : previously computed EFM for which we want to improve the accuracy (before removing the zero entries)
%       - A_lin : matrix which implements the linear constraints A_lin*x >= b_lin where x is the variable vector
%       - b_lin : vector which implements the linear constraints A_lin*x >= b_lin where x is the variable vector
%       - A_eq : matrix which implements the linear constraints A_eq*x = b_eq where x is the variable vector
%       - f : vector which implements the linear objective function f'*x where x is the variable vector  
% 
%   @outputs:  
%       - e_opt_new : new EFM matrix for which the accuracy has possibly been improved

  options_lin = optimoptions('linprog','Display','none');
  options_lin.Preprocess = 'none';
  ne = length(e_opt); e_opt_new = zeros(ne,1);
  ind_e_nonzeros = nonzeros(diag(abs(e_opt) > 1e-9)*(1:ne)');
  A_lin_red = A_lin(:,ind_e_nonzeros);
  A_eq_red = A_eq(:,ind_e_nonzeros);
  f_red = f(ind_e_nonzeros);
  e_opt_red = linprog(f_red,-A_lin_red,-b_lin,A_eq_red,b_eq,[],[],options_lin);
  e_opt_new(ind_e_nonzeros) = e_opt_red;

end