%% Build a matrix of double component terms, for efficiency purposes
%   author        : Mirko Pasquini 
%
%  In evaluating sensitivities (with respect to concentrations or
%  parameters) it happens that many of the same double components are
%  evaluated for several terms in the sensitivity matrices. This function
%  pre-calculate these terms and returns them stored in a matrix, to be
%  accessed to increase efficiency of the functions requiring them. 
%
%  NOTE: This efficiency gain might be negligible, tests to be done.
%
%  HDC_matrix = BuildHDCMatrix(thetaMatrix,c)
%
%  @inputs:
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.
%       c : vector of extracellular metabolite concentrations (> 0). 
%           These corresponds to coupled metabolites 
%           (Please refer to the "Beginner's guide" for a definition of 
%           coupled metabolite")
%
%   @outputs:
%       HDC_matrix : a matrix of double component Monod terms. 
%           If ka_ij is the activation coefficient for effect of the j-th 
%           metabolite on the i-th macroreaction rate, and ki_ij is the 
%           inhibition coefficient for the effect of the j-th metabolite on 
%           the i-th macroreaction rate, then the element (i,j) of the
%           matric HDC_matrix is hdc(c(j),ka_ij,ki_ij).

function HDC_matrix = BuildHDCMatrix(theta_matrix,c)
    HDC_matrix = zeros(size(theta_matrix,1),(size(theta_matrix,2)-1)/2);
    for i = 1 : size(theta_matrix,1)
        for j = 1 : (size(theta_matrix,2)-1)/2
            HDC_matrix(i,j) = hdc(c(j),theta_matrix(i,2*j),theta_matrix(i,2*j+1));
            if isinf(HDC_matrix(i,j))
                disp('Element ('+num2str(i)+','+num2str(j)+') of the HDC_matrix is inf.');
            end
        end
    end
end