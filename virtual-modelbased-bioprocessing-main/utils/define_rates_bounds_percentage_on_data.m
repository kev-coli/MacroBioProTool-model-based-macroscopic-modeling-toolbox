%% Define rates bounds based on a percentage of the data
%  author   :   Mirko Pasquini 
%
%   bounds = DEFINE_RATES_BOUNDS_PERCENTAGE_ON_DATA(qext_data,p)
%
%   Given a matrix of rates data-points it defines a "percentage bound"
%   around the data, i.e. a box containing all the data points, with the
%   lower bounds smaller (of a given percentage) than the minimum rate in the
%   data, and the upper bound larger (of a given percentage) than the maximum
%   rate in the data. The rationale behind this is that in a model-based
%   optimization if these bounds are satisfied, we are sure to be "close
%   enough" to the data on which the model is constructed, for higher trust
%   on the optimization results.
%
%   @inputs:
%       qext_data : matrix of rates measurements data, where each column
%           corresponds to a different tested condition and each row
%           corresponds to a particular extracellular metabolite (check the
%           stoichiometric matrix of the particular model for metabolite 
%           order);
%       p : percentage of variation for the defined bound. E.g. if p =
%           0.1, then lower bound will be 1-p = 0.9 of the minimum rates, 
%           while the upper bound will be 1+p = 1.1 of the maximum rates.
%
%   @outputs:
%       bounds : upper and lower bounds for the rates, based on the data
%           and the percentages specified.


function bounds = define_rates_bounds_percentage_on_data(qext_data,p)
    bounds = zeros(size(qext_data,1),2);
    min_qext = min(qext_data,[],2);
    max_qext = max(qext_data,[],2);
    for i = 1 : size(bounds,1)
        if min_qext(i) <= 0
            bounds(i,1) = min_qext(i)*(1+p);
        else
            bounds(i,1) = min_qext(i)/(1+p);
        end
        if max_qext(i) <= 0
            bounds(i,2) = max_qext(i)/(1+p);
        else
            bounds(i,2) = max_qext(i)*(1+p);
        end
    end
end