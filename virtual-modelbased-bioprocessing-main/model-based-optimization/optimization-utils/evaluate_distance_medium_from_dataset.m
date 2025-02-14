%% Evaluate distance of a medium from a given set of media data
%   author  :   Mirko Pasquini 
%
%   The function evaluates the distance of a medium from a given dataset.
%   In particular the distance will be the shortest distance between the
%   medium and any other medium in the dataset. The user has a choice to
%   normalize this distance (based on the medium bounds) and choose the
%   type of norm used for the distance evaluation (1-norm or 2-norm).
%
%   [distance, index_nearest_medium] 
%   = EVALUATE_DISTANCE_MEDIUM_FROM_DATASET(medium,...
%                                           media_data,...
%                                           norm_type,...
%                                           normalization_flag,...
%                                           medium_bounds)
%
%   @inputs:
%       medium : vector containing the components of the medium relative to
%           decision variables of the optimization problem
%       media_data : matrix in which each column corresponds to a medium in 
%           the data. Only the components relative to decision variables of
%           the optimization problem are considered;
%       norm_type : type of norm considered to evaluate the distance. The
%           available options are the same as for the function 'norm' (see
%           "help norm" for further info), however common choices are '1' 
%           for the 1-norm and '2' for the 2-norm;
%       normalization_flag : if the flag is set to 1, all media will be
%           normalized so that their components are between 0 and 1. This
%           is done to reduce the effect of having different interval
%           lengths for different metabolites;
%       medium_bounds : a matrix of dimension (nu x 2), where nu is the 
%           number of metabolites indexed by 
%           index_data.decision_metabolite_indices with the first column
%           being the lower bounds of the decision variables medium, and 
%           the second column being their upper bounds.
%
%   @outputs:
%       distance : distance of the medium from the available media in the
%           dataset media_data. This is evaluated as the minimum distance
%           between the medium and every other medium in the data;
%       index_nearest_medium : index of the column in media_data
%           corresponding to the vector nearest to medium (i.e. the nearest
%           medium available in the data)
%




function [distance, index_nearest_medium] = evaluate_distance_medium_from_dataset(medium,...
                                                                                media_data,...
                                                                                norm_type,...
                                                                                normalization_flag,...
                                                                                medium_bounds)
    n_media_in_data = size(media_data,2);
    if normalization_flag == 1
        % all media will be normalized so that all the components are
        % between 0 and 1
        medium_range = medium_bounds(:,2) - medium_bounds(:,1);
        medium = (medium - medium_bounds(:,1))./medium_range; 
        for n = 1 : n_media_in_data
            media_data(:,n) = (media_data(:,n)-medium_bounds(:,1))./medium_range;
        end
    end
    
    distance_vector = zeros(n_media_in_data,1);
    for n = 1 : n_media_in_data
        distance_vector(n) = norm(medium-media_data(:,n),norm_type);
    end

    [distance_vector_sorted, media_indexes_sorted_by_distance] = sort(distance_vector);
    distance = distance_vector_sorted(1);
    index_nearest_medium = media_indexes_sorted_by_distance(1);
end
