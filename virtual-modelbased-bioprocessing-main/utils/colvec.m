%% Return the column vector version of the vector v
%   author        : Mirko Pasquini
%
%   The function takes a vector v and returns it in column form. If v is
%   already a column the resulting vector is unaltered.
%   
%   vc = colvec(v)
%
%   @inputs:
%       v : any row or column vector (if a matrix is entered as input,
%           it displays an error message)
%   @outputs:
%       vc : the vector v in column form (i.e. v' if v is a row, v if v
%           is already a column)

function vc = colvec(v)
    if size(v,1) == 1
        vc = v';
    elseif size(v,2) == 1
        vc = v;
    elseif size(v,1) ~= 1  && size(v,2) ~= 1
        disp('The variable is not a vector.');
    end
end