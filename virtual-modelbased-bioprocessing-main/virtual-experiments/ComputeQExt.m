%% Evaluate extracellular metabolite rates given concentration and model
%   author        : Mirko Pasquini 
%
%   The function takes as input a concentration vector, a matrix of
%   kinetics parameters and a macroreaction stoichiometric matrix, and 
%   return the corresponding extracellular metabolite rates, as well as the
%   macroreaction rates associated to the model EFMs.
%
%   [q,w] = ComputeQExt(c, theta, Amac)
%
%   @inputs
%       c : vector of extracellular metabolite concentrations (>= 0). 
%           These corresponds to coupled metabolites (please refer to the
%           "Beginner's guide" for a definition of coupled metabolite")
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.
%       Amac : macroreaction stoichiometric matrix, obtained as Aext*E, 
%           where Aext is the extracellular stoichiometric matrix and E is 
%           the matrix of EFMs.
%
%   @outputs
%       q : predicted extracellular metabolites rates given the model 
%           described by theta and Amac
%       w : macroreaction rates given the model's EFMs

function [q,w] = ComputeQExt(c,theta,Amac)
    m = (size(theta,2)-1)/2;
    w = zeros(size(Amac,2),1);
    for i = 1 : length(w)
        w(i) = theta(i,1);
        for j = 1 : length(c)
            w(i) = w(i)*hdc(c(j),theta(i,j*2),theta(i,j*2+1));
        end
    end
    q = Amac*w;
end