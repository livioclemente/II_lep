clear; clc; close all;

addpath('/Users/livioclemente/Documents/MATLAB/bluewhitered/');
addpath('/Users/livioclemente/Documents/MATLAB/eeglab2023.0/');
addpath('/Users/livioclemente/Documents/MATLAB/EKG_EEG-master/');
addpath('/Users/livioclemente/Documents/MATLAB/DrosteEffect_BrewerMap_3/');

time_start = -0.1; % start time in seconds
time_end = 0.6; % end time in seconds
srate = 256;
time_epoch = time_start:1/srate:time_end-1/srate;
n_samples = round((time_end - time_start) * srate);
time_info = linspace(time_start * 1000, time_end * 1000 - 1 / srate * 1000, n_samples);

% Tracciati per il ginocchio
data_path_knee_ct = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/MI_LEP/conditions/pth/piede/';
sub_knee_ct = dir([data_path_knee_ct '*.set']); % List of subject files for knee
nsub_knee_ct = length(sub_knee_ct); % numero soggetti

% Preallocazione dei cell arrays
DATA_knee_ct = cell(nsub_knee_ct, 1);
x_values_knee_ct = cell(nsub_knee_ct, 1);
LEP_knee_ct = cell(nsub_knee_ct, 1);

for isub = 1:nsub_knee_ct % Carica i dati completi
    EEG_knee_ct = pop_loadset(sub_knee_ct(isub).name, data_path_knee_ct);
    DATA_knee_ct{isub} = EEG_knee_ct.data;

    [electrodes, Nt, trials_knee] = size(DATA_knee_ct{isub}); % Controlla la dimensione dei dati caricati
    
    % Calcolare l'indice dei campioni corrispondenti a -0.1 e 0.6 secondi
    sample_start = 1; % -0.1 sec in samples
    sample_end = round((time_end + abs(time_start)) * srate); % 0.6 sec in samples

    % Controlla che sample_end non superi il numero di campioni disponibili
    if sample_end > Nt
        sample_end = Nt;  % Limita sample_end al numero massimo di campioni
    end
    
    % LEP
    LEP_knee_ct{isub} = mean(DATA_knee_ct{isub}(:, sample_start:sample_end, :), 3);
end

LEP_avg_knee_ct = mean(cat(3, LEP_knee_ct{:}), 3);

idx_n2 = find(time_info >= 190 & time_info <= 250);
idx_p2 = find(time_info >= 310 & time_info <= 395);

% Preallocazione delle variabili per memorizzare i valori di ampiezza e latenza
n2_amp_values = zeros(nsub_knee_ct, 1);
n2_latency_values = zeros(nsub_knee_ct, 1);
p2_amp_values = zeros(nsub_knee_ct, 1);
p2_latency_values = zeros(nsub_knee_ct, 1);

for isub = 1:nsub_knee_ct
    % Estrai il LEP medio per il soggetto corrente
    lep_signal = LEP_knee_ct{isub}(11, :);

    % Calcola ampiezza e latenza per N2
    [n2_amp, n2_idx] = min(lep_signal(idx_n2)); % Supponendo che N2 sia un picco negativo
    n2_latency = time_info(idx_n2(n2_idx));
    n2_amp_values(isub) = n2_amp;
    n2_latency_values(isub) = n2_latency;

    % Calcola ampiezza e latenza per P2
    [p2_amp, p2_idx] = max(lep_signal(idx_p2)); % Supponendo che P2 sia un picco positivo
    p2_latency = time_info(idx_p2(p2_idx));
    p2_amp_values(isub) = p2_amp;
    p2_latency_values(isub) = p2_latency;
end

mean_n2_amp = mean(n2_amp_values);
mean_n2_latency = mean(n2_latency_values);
mean_p2_amp = mean(p2_amp_values);
mean_p2_latency = mean(p2_latency_values);

% Verifica che i valori non siano NaN o 0
disp(['Mean N2 Latency: ', num2str(mean_n2_latency)]);
disp(['Mean P2 Latency: ', num2str(mean_p2_latency)]);
disp(['Mean N2 Amplitude: ', num2str(mean_n2_amp)]);
disp(['Mean P2 Amplitude: ', num2str(mean_p2_amp)]);


