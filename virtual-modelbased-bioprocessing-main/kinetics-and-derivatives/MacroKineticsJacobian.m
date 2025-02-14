%% Jacobian of macroreaction kinetics
%   author        : Mirko Pasquini 
%
%  ======================================================================
%  There might be some formal definition problems if a concentration is
%  exactly 0. For this reason it is always a good idea to preprocess the
%  concentrations to be strictly positive e.g. 1e-4.
%  ======================================================================
%
%  The function evaluates the Jacobian of the macroreaction kinetic
%  functions. The function takes as input a parameters matrix and a vector
%  of concentrations, in which the Jacobian needs to be evaluated, and it
%  returns the Jacobian. 
%
%  J_w = MacroKineticsJacobian(theta_matrix,c,macro_of_interest)
%
%  @inputs:
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.
%       c : vector of extracellular metabolite concentrations (> 0). 
%           These corresponds to coupled metabolites (Please refer to the 
%           "Beginner's guide"for a definition of coupled metabolite")
%       macro_of_interest : a vector of 0s and 1s of the same size of the
%           macrorates vector. If macro_of_interest(i) = 0,
%           then we are not interested in the gradient of macrorate w_i, so
%           we can avoid the computation. If this argument is not passed,
%           the complete Jacobian is evaluated.
%
%
%  @outputs:
%       J_w : Jacobian of the macroreaction kinetics with respect to the 
%           concentrations c. In particular the element (i,j) of such 
%           matrix would be dw_i/dc_j


function J_w = MacroKineticsJacobian(theta_matrix,c,macro_of_interest)
    if nargin < 3
        macro_of_interest = ones(size(theta_matrix,1));
    end
    Ka = [];
    Ki = [];
    WM = theta_matrix(:,1);
    
    for i = 1 : length(c)
        Ka = [Ka,theta_matrix(:,i*2)];      
        Ki = [Ki,theta_matrix(:,i*2+1)];   
    end
    
    r = size(theta_matrix,1);
    m = length(c);
    J_w = Jacobian_of_kinetic_macro(Ka,Ki,WM,r,m,c,macro_of_interest);  
end

function J = Jacobian_of_kinetic_macro(Ka,Ki,WMAX,r,m,y,macro_of_interest)
    J = zeros(r,m);
    for i = 1  : r
        if macro_of_interest(i) ~= 0
            for j = 1 : m
                J_ij = WMAX(i);
                J_ij = J_ij*(partial_kin_yj_wi(Ka(i,j),Ki(i,j),y(j)));
                Ka_minus_j = Ka;
                Ka_minus_j(:,j) = [];
                Ki_minus_j = Ki;
                Ki_minus_j(:,j) = [];
                y_minus_j = y;
                y_minus_j(j) = [];
                J_ij = J_ij*(double_comp_multiplication_minus_yj(Ka_minus_j(i,:),Ki_minus_j(i,:),y_minus_j));
                J(i,j) = J_ij;
            end
        end
    end
end

function dhdj = partial_kin_yj_wi(Ka_i_j,Ki_i_j,y_j)
    if Ka_i_j == 0 && Ki_i_j == 0
        dhdj = 0;
    elseif Ki_i_j == 0 % act
        dhdj = (Ka_i_j)/(y_j + Ka_i_j)^2;
    elseif Ka_i_j == 0 % inh
        dhdj = -(Ki_i_j)/(1+y_j*Ki_i_j)^2;
    else
        dhdj = (Ka_i_j - y_j^2*Ki_i_j)/((y_j+Ka_i_j)^2*(1+y_j*Ki_i_j)^2);
    end
end

function dcmj = double_comp_multiplication_minus_yj(Ka_minus_j_i,Ki_minus_j_i,y_minus_j)
    dcmj = 1;
    for k = 1 : length(y_minus_j)
        factor = hdc(y_minus_j(k), Ka_minus_j_i(k), Ki_minus_j_i(k)); 
        dcmj = dcmj*factor; 
    end
end