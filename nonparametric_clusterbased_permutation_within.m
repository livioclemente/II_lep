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

% Carica i valori VAS da un file Excel
vas_file = [data_path_knee_ct 'vas.xlsx'];
vas_values = readmatrix(vas_file);
if length(vas_values) ~= nsub_knee_ct
    error('Number of VAS values does not match the number of subjects');
end

%% Calcola le informazioni temporali per soggetto

% Preallocazione dei cell arrays
DATA_knee_ct = cell(nsub_knee_ct, 1);
x_values_knee_ct = cell(nsub_knee_ct, 1);
vas_values_knee_ct = cell(nsub_knee_ct, 1);
MI_knee_ct = cell(nsub_knee_ct, 1);
II_knee_ct = cell(nsub_knee_ct, 1);
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

    % Taglia i dati solo prima di calcolare MI e II
    x = squeeze(DATA_knee_ct{isub}(11, sample_start:sample_end, :))';
    vas = repmat(vas_values(isub), size(x, 1), 1); % Ripeti il valore VAS per ogni trial per il soggetto corrente
    x_knee_ct = copnorm(x);
    vas_knee_ct = copnorm(vas);

    x_values_knee_ct{isub} = x_knee_ct;
    vas_values_knee_ct{isub} = vas_knee_ct;

    MI = zeros(1, size(x, 2));
    for ti = 1:size(x, 2)
        MI(ti) = mi_gg(x_knee_ct(:, ti), vas_knee_ct(:, 1), true, true);
    end
    
    MI_knee_ct{isub} = MI;

    noise = .00000005 * randn(size(x_knee_ct, 1), 1);
    II = zeros(size(x, 2), size(x, 2));
    for t1 = 1:size(x, 2)
        for t2 = (t1 + 1):size(x, 2)
            JMI = mi_gg([x_knee_ct(:, t1) x_knee_ct(:, t2) + noise], vas_knee_ct(:, 1), true, true);
            II(t1, t2) = JMI - MI(t1) - MI(t2);
        end
    end
    II = II + II';
    II_knee_ct{isub} = II;
    
    % LEP
    LEP_knee_ct{isub} = mean(DATA_knee_ct{isub}(:, sample_start:sample_end, :), 3);
end

%% Global

MI_avg_knee_ct = mean(cat(1, MI_knee_ct{:}), 1);
II_avg_knee_ct = mean(cat(3, II_knee_ct{:}), 3);
LEP_avg_knee_ct = mean(cat(3, LEP_knee_ct{:}), 3);

%% Calcolo Aggregato e Permutazioni

% Preallocazione per aggregazione
n_trials_total = sum(cellfun(@(x) size(x, 3), DATA_knee_ct)); % Numero totale di trials
n_timepoints = sample_end - sample_start + 1; % Numero di timepoints per trial
all_x_knee_ct = zeros(n_trials_total, n_timepoints); % Preallocazione per tutti i trials e timepoints
all_vas_knee_ct = zeros(n_trials_total, 1); % Preallocazione per VAS normalizzati
current_idx = 1;


% Caricamento e aggregazione dei dati
for isub = 1:nsub_knee_ct
    EEG_knee_ct = pop_loadset(sub_knee_ct(isub).name, data_path_knee_ct);
    
    % Taglia i dati solo prima di calcolare MI e II
    x = squeeze(EEG_knee_ct.data(11, sample_start:sample_end, :))'; % Estrazione dati per canale specifico
    vas = repmat(vas_values(isub), size(x, 1), 1); % Ripeti il valore VAS per ogni trial
    
    % Normalizzazione e aggregazione
    n_trials = size(x, 1); % Numero di trials per il soggetto corrente
    all_x_knee_ct(current_idx:current_idx + n_trials - 1, :) = copnorm(x); % Inserisce i dati normalizzati nel posto giusto
    all_vas_knee_ct(current_idx:current_idx + n_trials - 1) = copnorm(vas); % Inserisce VAS normalizzati
    current_idx = current_idx + n_trials; % Aggiorna l'indice corrente
end

% Calcolo della MI aggregata
MI_aggregated = zeros(1, size(all_x_knee_ct, 2));
for ti = 1:size(all_x_knee_ct, 2)
    MI_aggregated(ti) = mi_gg(all_x_knee_ct(:, ti), all_vas_knee_ct, true, true);
end

nperm = 2000; % number of permutations
preclust_pval = 0.05;
clust_pval = 0.05;

% Calcolo dell'II aggregata
II_aggregated = zeros(size(all_x_knee_ct, 2), size(all_x_knee_ct, 2));
noise = .00000005 * randn(size(all_x_knee_ct, 1), 1);
h0 = zeros(size(all_x_knee_ct, 2), size(all_x_knee_ct, 2));
hp = zeros(size(all_x_knee_ct, 2), size(all_x_knee_ct, 2), nperm);
clust_max = zeros(nperm,1);

for t1 = 1:size(all_x_knee_ct, 2)
    for t2 = (t1 + 1):size(all_x_knee_ct, 2)
        JMI_aggregated = mi_gg([all_x_knee_ct(:, t1) all_x_knee_ct(:, t2) + noise], all_vas_knee_ct, true, true);
        II_aggregated(t1, t2) = JMI_aggregated - MI_aggregated(t1) - MI_aggregated(t2);
        
        h0(t1, t2) = II_aggregated(t1, t2);
        
        for iperm = 1:nperm
            vas_perm = all_vas_knee_ct(randperm(length(all_vas_knee_ct))); % permute VAS values
            JMI_perm = mi_gg([all_x_knee_ct(:, t1) all_x_knee_ct(:, t2) + noise], vas_perm, true, true);
            II_perm = JMI_perm - mi_gg(all_x_knee_ct(:, t1), vas_perm, true, true) - mi_gg(all_x_knee_ct(:, t2), vas_perm, true, true);
            hp(t1, t2, iperm) = II_perm;
        end
    end
end

% fill in the other half of the matrices
h0 = h0 + h0';
for iperm = 1:nperm 
    hp(:,:,iperm) = hp(:,:,iperm) + hp(:,:,iperm)'; 
end 

%% collect the largest suprathreshold clusters
for iperm = 1:nperm
    perms = true(1,nperm);
    perms(iperm) = 0;
    zvals = squeeze((hp(:,:,iperm) - mean(hp(:,:,perms),3)) ./ std(hp(:,:,perms),[],3)); % see Cohen's book chapter 33 equation 33.1
    zvals(abs(zvals) < norminv(1 - preclust_pval)) = 0; 
    
    clust_info = bwconncomp(zvals);
    clust_max(iperm) = max([0 cellfun(@numel,clust_info.PixelIdxList)]);
end

%% identify significant clusters in real data
zmap = squeeze((h0 - mean(hp,3)) ./ std(hp,[],3));
zmap(abs(zmap) < norminv(1 - preclust_pval)) = 0;

clust_info = bwconncomp(zmap);
clust_size = cellfun(@numel,clust_info.PixelIdxList);
clust_th = prctile(clust_max,100-clust_pval*100);
clust_rem = find(clust_size < clust_th);
for i=1:length(clust_rem)
    zmap(clust_info.PixelIdxList{clust_rem(i)})=0;
end
zmap(isnan(zmap)) = 0;
zmap = logical(zmap);

%% plot
figure;
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
time_plot = time_info*1000;
imagesc(time_plot,time_plot,h0(:,:))
hold on
contour(time_plot,time_plot,zmap,1,'linecolor','k','LineWidth',1)
colormap(brewermap([],'*RdBu'));
colorbar
axis square
lim = max(abs(min(min(h0(:,:)))),max(max(h0(:,:))));
clim([-lim,lim]) % center zero

ax2 = subplot(5, 5, [22.2 23 24 24.35]);
% Traccia il LEP medio per il canale specifico (canale 11)
plot(time_plot, LEP_avg_knee_ct(11, :), 'k', 'LineWidth', 2); % Modificato il colore della linea in rosso per maggiore visibilità
axis tight;
xlim([min(time_plot), max(time_plot)]);
xlabel('Time (ms)'),ylabel('Time (ms)')
box off;


%% plot con maschera
figure;
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
time_plot = time_info * 1000;

% Applica la maschera: imposta a 0 tutto ciò che non è significativo
h0_masked = h0;
h0_masked(~zmap) = 0; % Imposta a 0 i valori fuori dai cluster significativi

% Visualizza la matrice modificata
imagesc(time_plot, time_plot, h0_masked);
hold on;
contour(time_plot, time_plot, zmap, 1, 'linecolor', 'k', 'LineWidth', 1);
colormap(brewermap([], '*RdBu'));
colorbar;
axis square;
lim = max(abs(min(min(h0(:,:)))), max(max(h0(:,:))));
caxis([-lim, lim]); % center zero

ax2 = subplot(5, 5, [22.2 23 24 24.35]);
% Traccia il LEP medio per il canale specifico (canale 11)
plot(time_plot, LEP_avg_knee_ct(11, :), 'k', 'LineWidth', 2); % Modificato il colore della linea in nero per maggiore visibilità
axis tight;
xlim([min(time_plot), max(time_plot)]);
xlabel('Time (ms)');
ylabel('Time (ms)');
box off;

clear

%% clusters info - Creazione della struttura "clusters" per memorizzare i dettagli dei cluster significativi

% Ricalcola la mappa dei cluster dopo aver applicato la soglia di significatività
clust_info_signif = bwconncomp(zmap);

% Numero di cluster significativi
num_signif_clusters = clust_info_signif.NumObjects;

% Calcola la somma dei valori t per i cluster significativi
clust_sum_signif = cellfun(@(x) sum(h0(x)), clust_info_signif.PixelIdxList);

% Calcola la dimensione di ciascun cluster significativo
clust_size_signif = cellfun(@numel, clust_info_signif.PixelIdxList);

clusters = struct(); % Definisci la struttura iniziale vuota

for i = 1:num_signif_clusters
    % Indici dei punti che compongono il cluster i
    cluster_indices = clust_info_signif.PixelIdxList{i};

    % Converti gli indici lineari in indici di matrice per ottenere i punti temporali
    [t1_indices, t2_indices] = ind2sub(size(zmap), cluster_indices);

    % Ottieni i tempi corrispondenti
    t1_times = time_info(t1_indices);
    t2_times = time_info(t2_indices);

    % Peso del cluster (somma dei valori t)
    cluster_weight = sum(h0(cluster_indices)); % Usa h0 per ottenere i valori t originali

    % Dimensione del cluster (numero di punti che compongono il cluster)
    cluster_size = numel(cluster_indices);

    % Memorizza le informazioni nella struttura
    clusters(i).number = i;
    clusters(i).weight = cluster_weight;
    clusters(i).size = cluster_size;
    clusters(i).t1_min = min(t1_times);
    clusters(i).t1_max = max(t1_times);
    clusters(i).t2_min = min(t2_times);
    clusters(i).t2_max = max(t2_times);
    clusters(i).t1_times = t1_times;
    clusters(i).t2_times = t2_times;
    
    % (Opzionale) Memorizza i tempi completi
    clusters(i).t1_times = t1_times;
    clusters(i).t2_times = t2_times;
    
    % Visualizza le informazioni (opzionale per debugging)
    fprintf('Cluster %d:\n', clusters(i).number);
    fprintf(' - Peso (somma dei valori t): %.2f\n', clusters(i).weight);
    fprintf(' - Dimensione: %d punti\n', clusters(i).size);
    fprintf(' - Intervalli temporali coinvolti:\n');
    fprintf('   * Time 1: %.2f ms to %.2f ms\n', clusters(i).t1_min, clusters(i).t1_max);
    fprintf('   * Time 2: %.2f ms to %.2f ms\n\n', clusters(i).t2_min, clusters(i).t2_max);
end

% Calcolo dei valori di p e t per ogni confronto temporale, se disponibili
p_values = nan(size(h0));
t_values = nan(size(h0));

% Calcolo e salvataggio dei valori p e t
for t1 = 1:size(h0, 1)
    for t2 = (t1 + 1):size(h0, 2)
        % Estrai le variabili originali per calcolare il test t
        data_ct = all_x_knee_ct(:, t1); % Dati di un gruppo
        data_pth = all_x_knee_ct(:, t2); % Dati dell'altro gruppo
        [~, p, ~, stats] = ttest2(data_ct, data_pth);
        p_values(t1, t2) = p;
        t_values(t1, t2) = stats.tstat;
    end
end

% Completa le matrici p_values e t_values
p_values = p_values + p_values';
t_values = t_values + t_values';

% Prealloca variabili per le statistiche dei cluster
clust_mean_t = zeros(clust_info_signif.NumObjects, 1);
clust_std_t = zeros(clust_info_signif.NumObjects, 1);

% Calcola la media, la deviazione standard, e la distribuzione spaziale per ciascun cluster
for i = 1:clust_info_signif.NumObjects
    cluster_indices = clust_info_signif.PixelIdxList{i};
    clust_values = h0(cluster_indices); % Prendi i valori t/z del cluster
    clust_mean_t(i) = mean(clust_values);
    clust_std_t(i) = std(clust_values);
end
%% Salvataggio delle Variabili
% Aggiungi la struttura dei cluster al salvataggio delle variabili principali
variables_to_save = {'h0', 'zmap', 'time_info', 'lim', 'clust_info_signif', 'clust_sum_signif', 'clust_size_signif', 'clusters', 'LEP_avg_knee_ct', 'hp', 'p_values', 't_values', 'clust_mean_t', 'clust_std_t'};

% Salva tutte le variabili in un file .mat
save('2kperm_within_pth_piede.mat', variables_to_save{:});
