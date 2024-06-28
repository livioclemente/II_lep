clear; clc;

addpath('/Users/livioclemente/Documents/MATLAB/eeglab2023.0/');
addpath('/Users/livioclemente/Documents/MATLAB/EKG_EEG-master/');
addpath('/Users/livioclemente/Documents/MATLAB/bluewhitered/');

% Definition of groups and districts
groups = {'ct', 'dl_pth', 'normali', 'pth'};
districts = {'ginocchio', 'mano', 'piede'};
base_path = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/MI_LEP/conditions/';

time_start = -0.1; % start time in seconds
time_end = 1; % end time in seconds
srate = 256;
time_epoch = time_start:1/srate:time_end-1/srate;
n_samples = round((time_end - time_start) * srate);
time_info = linspace(time_start * 1000, time_end * 1000 - 1 / srate * 1000, n_samples);
ntime_info = length(time_info);

%% Calculate temporal information per group per district per subject

for g = 1:length(groups)
    group = groups{g};
    fprintf('Processing group: %s\n', group);

    for d = 1:length(districts)
        district = districts{d};
        fprintf('Processing district: %s in group: %s\n', district, group);
        data_path = [base_path group '/' district '/'];
        var_prefix = sprintf('%s_%s', group, district); % Prefix for variable names
        eval(sprintf('%s_sub_files = dir([data_path ''*.set'']);', var_prefix));
        eval(sprintf('%s_nsub = length(%s_sub_files);', var_prefix, var_prefix));

        % Preallocation of cell arrays
        eval(sprintf('%s_DATA = cell(%s_nsub, 1);', var_prefix, var_prefix));
        eval(sprintf('%s_x = cell(%s_nsub, 1);', var_prefix, var_prefix));
        eval(sprintf('%s_vas = cell(%s_nsub, 1);', var_prefix, var_prefix));  % Preallocation for vas
        eval(sprintf('%s_MI = cell(%s_nsub, 1);', var_prefix, var_prefix));
        eval(sprintf('%s_II = cell(%s_nsub, 1);', var_prefix, var_prefix));
        eval(sprintf('%s_LEP = cell(%s_nsub, 1);', var_prefix, var_prefix));

        % Preallocazione per aggregazione
        eval(sprintf('%s_aggregated_x = [];', var_prefix));
        eval(sprintf('%s_aggregated_vas = [];', var_prefix));

        % Load vas data from Excel file
        vas_file = [data_path 'vas.xlsx'];
        vas_values = readmatrix(vas_file);
        if length(vas_values) ~= eval(sprintf('%s_nsub', var_prefix))
            error('Number of VAS values does not match the number of subjects');
        end

        % Load data    
        for isub = 1:eval(sprintf('%s_nsub', var_prefix))
            EEG = pop_loadset(eval(sprintf('%s_sub_files(isub).name', var_prefix)), data_path);
            eval(sprintf('%s_DATA{isub} = EEG.data;', var_prefix));

            % Extraction and normalization of data
            [electrodes, Nt, trials] = size(eval(sprintf('%s_DATA{isub}', var_prefix))); 
            x = squeeze(eval(sprintf('%s_DATA{isub}(11, :, :)', var_prefix)))';
            vas = repmat(vas_values(isub), trials, 1);  % Repeat vas_value for each trial for the current subject
            x_values = copnorm(x);
            vas_values_norm = copnorm(vas);

            eval(sprintf('%s_x{isub} = x_values;', var_prefix));
            eval(sprintf('%s_vas{isub} = vas_values_norm;', var_prefix));

            % Aggregate data
            eval(sprintf('%s_aggregated_x = [%s_aggregated_x; x_values];', var_prefix, var_prefix));
            eval(sprintf('%s_aggregated_vas = [%s_aggregated_vas; vas_values_norm];', var_prefix, var_prefix));

            % MI
            MI = zeros(1, Nt);
            for ti = 1:Nt
                MI(ti) = mi_gg(x_values(:, ti), vas_values_norm(:, 1), true, true);
            end
            eval(sprintf('%s_MI{isub} = MI;', var_prefix));

            % Calculate JMI and II
            noise = .00000005 * randn(size(x_values, 1), 1);
            II = zeros(Nt, Nt);
            for t1 = 1:Nt
                for t2 = (t1+1):Nt
                    JMI = mi_gg([x_values(:, t1) x_values(:, t2) + noise], vas_values_norm(:, 1), true, true);
                    II(t1, t2) = JMI - MI(t1) - MI(t2);
                end
            end
            II = II + II';
            eval(sprintf('%s_II{isub} = II;', var_prefix));

            % MI per canale
            MI_all_channels = zeros(Nt, electrodes);
            for ci = 1:electrodes
                x_channel = squeeze(eval(sprintf('%s_DATA{isub}(ci, :, :)', var_prefix)))';
                x_values_channel = copnorm(x_channel);
                for ti = 1:Nt
                    MI_all_channels(ti, ci) = mi_gg(x_values_channel(:, ti), vas_values_norm(:, 1), true, true);
                end
            end
            eval(sprintf('%s_MI_all_channels{isub} = MI_all_channels;', var_prefix));

            % LEP
            eval(sprintf('%s_LEP{isub} = mean(%s_DATA{isub}, 3);', var_prefix, var_prefix));
        end

        % Aggregate MI and II calculation
        eval(sprintf('%s_MI_aggregated = zeros(1, size(%s_aggregated_x, 2));', var_prefix, var_prefix));
        for ti = 1:size(eval(sprintf('%s_aggregated_x', var_prefix)), 2)
            eval(sprintf('%s_MI_aggregated(ti) = mi_gg(%s_aggregated_x(:, ti), %s_aggregated_vas, true, true);', var_prefix, var_prefix, var_prefix));
        end

        eval(sprintf('%s_II_aggregated = zeros(size(%s_aggregated_x, 2), size(%s_aggregated_x, 2));', var_prefix, var_prefix, var_prefix));
        noise = .00000005 * randn(size(eval(sprintf('%s_aggregated_x', var_prefix)), 1), 1);
        for t1 = 1:size(eval(sprintf('%s_aggregated_x', var_prefix)), 2)
            for t2 = (t1+1):size(eval(sprintf('%s_aggregated_x', var_prefix)), 2)
                eval(sprintf('JMI_aggregated = mi_gg([%s_aggregated_x(:, t1) %s_aggregated_x(:, t2) + noise], %s_aggregated_vas, true, true);', var_prefix, var_prefix, var_prefix));
                eval(sprintf('%s_II_aggregated(t1, t2) = JMI_aggregated - %s_MI_aggregated(t1) - %s_MI_aggregated(t2);', var_prefix, var_prefix, var_prefix));
            end
        end
        eval(sprintf('%s_II_aggregated = %s_II_aggregated + %s_II_aggregated'';', var_prefix, var_prefix, var_prefix));

        % Global values
        eval(sprintf('%s_MI_avg = mean(cat(1, %s_MI{:}), 1);', var_prefix, var_prefix));
        eval(sprintf('%s_II_avg = mean(cat(3, %s_II{:}), 3);', var_prefix, var_prefix));
        eval(sprintf('%s_LEP_avg = mean(cat(3, %s_LEP{:}), 3);', var_prefix, var_prefix));
        eval(sprintf('%s_MI_all_channels_avg = mean(cat(3, %s_MI_all_channels{:}), 3);', var_prefix, var_prefix));

        % Save
        save_file = sprintf('%s.mat', var_prefix);
        eval(sprintf('save(save_file, ''%s_MI'', ''%s_MI_avg'', ''%s_II'', ''%s_II_avg'', ''%s_LEP'', ''%s_LEP_avg'', ''Nt'', ''%s_sub_files'', ''time_info'', ''%s_vas'', ''%s_nsub'', ''%s_x'', ''%s_aggregated_x'', ''%s_aggregated_vas'', ''%s_MI_aggregated'', ''%s_II_aggregated'', ''%s_MI_all_channels'', ''%s_MI_all_channels_avg'');', var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix, var_prefix));
    end
end

