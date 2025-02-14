%% Test-unit function 2 : model predictions from model and data
%
%   The function takes as input a model and paths of datafiles and generate
%   predictions on steady-state concentrations and rates (using
%   mass-balance equation or using the model directly).


function [model_concentrations_predictions_mbe,...
            model_rates_predictions_mbe,...
            model_rates_predictions_direct,...
            qext_data,...
            cext_data,...
            cin_data] = get_model_predictions_testing(model, datafiles, index_data)

    % load qext data
    qext_path = datafiles.qext_data_file;
    [~,~,qext_path_extension] = fileparts(qext_path);
    if (strcmp(qext_path_extension,".xlsx") || strcmp(qext_path_extension,".xls"))
        qext_table = readtable(qext_path);
        qext_data = qext_table{:,3:end}'; % be careful of transpose
    elseif (strcmp(qext_path_extension,".mat"))
        qext_data = load(qext_path);
    else
        disp("Unknown extension of qext file.");
    end

    % load cext data
    cext_path = datafiles.cext_data_file;
    [~,~,cext_path_extension] = fileparts(cext_path);
    if (strcmp(cext_path_extension,".xlsx") || strcmp(cext_path_extension,".xls"))
        cext_table = readtable(cext_path);
        cext_data = cext_table{:,3:end}'; % be careful of transpose
    elseif (strcmp(cext_path_extension,".mat"))
        cext_data = load(cext_path);
    else
        disp("Unknown extension of cext file.");
    end

    cin_path = datafiles.cin_data_file;
    [~,~,cin_path_extension] = fileparts(cin_path);
    if (strcmp(cin_path_extension,".xlsx") || strcmp(cin_path_extension,".xls"))
        cin_table = readtable(cin_path);
        cin_data = cin_table{:,3:end}'; % be careful of transpose
    elseif (strcmp(cin_path_extension,".mat"))
        cin_data = load(cin_path);
    else
        disp("Unknown extension of cin file.");
    end

    model_concentrations_predictions_mbe = zeros(size(cext_data));
    model_rates_predictions_mbe = zeros(size(qext_data));
    model_rates_predictions_direct = zeros(size(qext_data));

    n_conditions = size(qext_data,2);
    if n_conditions ~= size(cext_data,2)
        fprintf("Number of conditions is different between concentrations and rates measurements datafile. Please check.")
    end


    for n = 1 : n_conditions
        [model_concentrations_predictions_mbe(:,n), model_rates_predictions_mbe(:,n)] =...
            run_virtual_experiment(cin_data(:,n), model, cext_data(:,n), index_data.coupled_met_indices, 'none', 1e-5, 500);
        model_rates_predictions_direct(:,n) = ComputeQExt(cext_data(:,n), model.theta_matrix, model.Amac);
    n
    end


end