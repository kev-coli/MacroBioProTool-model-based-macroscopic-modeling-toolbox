%% Return the row vector version of the vector v
%   author        : Mirko Pasquini 
%
%   The function takes a vector v and returns it in row form. If v is
%   already a row the resulting vector is unaltered.
%   
%   vc = rowvec(v)
%
%   @inputs:
%       v : any row or column vector (if a matrix is entered as input,
%           it displays an error message)
%   @outputs:
%       vr : the vector v in row form (i.e. v' if v is a column, v if v
%           is already a row)

function vr = rowvec(v)
    if size(v,2) == 1
        vr = v';
    elseif size(v,1) == 1
        vr = v;
    elseif size(v,1) ~= 1  && size(v,2) ~= 1
        disp('The variable is not a vector.');
    end
end