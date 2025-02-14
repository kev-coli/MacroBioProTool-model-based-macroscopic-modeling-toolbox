%% Evaluate the vector of sensitivities (w.r.t. kinetic parameters) of the 
%% i-th macrorate
%  author        : Mirko Pasquini 
%
%  The function evaluates the vector of sensitivities (w.r.t. kinetic parameters) 
%  of the i-th macrorate, where the parameters considered are the ones in
%  the i-th row of the parameters matrix (the sensitivity is zero with
%  respect to all the other macroreactions parameters, as they are
%  decoupled). The sensitivity depends on the particular extracellular
%  metabolites concentrations vector we are considering.
%
%  nabla_theta_wi = BUILD_ITH_MACRORATE_SENSITIVITY_VECTOR(HDC_matrix_i,theta_matrix_i,c)
%
%  @inputs:
%       HDC_matrix_i : i-th row of the matrix HDC_matrix, as returned by 
%           the function BuildHDCMatrix (please refer to that function for 
%           further details).
%       theta_matrix_i : i-th row of the kinetic parameters matrix. Please 
%           refer to the "Beginner's guide" for further info on the kinetic
%           parameters matrix.
%       c : vector of extracellular metabolite concentrations (> 0). These 
%           corresponds to coupled metabolites (please refer to the 
%           "Beginner's guide" for a definition of coupled metabolite").
%
%  @outpus:
%       nabla_theta_wi : vector of sensitivities of the i-th macroreaction 
%           with respect to its kinetic parameters. Its j-th element is
%           (dw_i/dtheta_matrix(i,j)).



function nabla_theta_wi = build_ith_macrorate_sensitivity_vector(HDC_matrix_i,theta_matrix_i,c)
    nabla_theta_wi = ones(1,length(theta_matrix_i));
    for j = 1 : length(theta_matrix_i)
        if j == 1
            for k = 1 : length(c)
                nabla_theta_wi(j) = nabla_theta_wi(j)*HDC_matrix_i(k);
            end
            %maximal rate sensitivity
        elseif mod(j,2) == 0
            %activation coeff sensitivity
            for k = 1 : length(c)
                if k ~= j/2
                    nabla_theta_wi(j) = nabla_theta_wi(j)*HDC_matrix_i(k);
                else
                    theta_act = theta_matrix_i(2*k);
                    theta_inh = theta_matrix_i(2*k+1);
                    nabla_theta_wi(j) = nabla_theta_wi(j)*hdc(c(k),0,theta_inh);
                    if c(k) == 0   
                        nabla_theta_wi(j) = 0;
                    else
                        nabla_theta_wi(j) = -nabla_theta_wi(j)*(c(k)/(theta_act+c(k))^2);
                    end
                end
            end
            nabla_theta_wi(j) = nabla_theta_wi(j)*theta_matrix_i(1);
        else
            %inhibition coeff sensitivity
            for k = 1 : length(c)
                if k ~= (j-1)/2
                    nabla_theta_wi(j) = nabla_theta_wi(j)*HDC_matrix_i(k);
                else
                    theta_act = theta_matrix_i(2*k);
                    theta_inh = theta_matrix_i(2*k+1);
                    nabla_theta_wi(j) = nabla_theta_wi(j)*hdc(c(k),theta_act,0);
                    if c(k) == 0
                        nabla_theta_wi(j) = 0;
                    else
                        nabla_theta_wi(j) = -nabla_theta_wi(j)*(c(k)/(1+theta_inh*c(k))^2);
                    end
                end
            end
            nabla_theta_wi(j) = nabla_theta_wi(j)*theta_matrix_i(1);
        end
        if isinf(nabla_theta_wi(j))
            disp('What is happened here?') % no entry should be infinity, this is just a check for debugging.
        end
    end
end