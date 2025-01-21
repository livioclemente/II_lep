clear; clc;

addpath('/Users/livioclemente/Documents/MATLAB/bluewhitered/');
addpath('/Users/livioclemente/Documents/MATLAB/eeglab2023.0/');
addpath('/Users/livioclemente/Documents/MATLAB/EKG_EEG-master/');
addpath('/Users/livioclemente/Documents/MATLAB/DrosteEffect_BrewerMap_3/');

% Carica i dati dei gruppi
load('ct_mano.mat');
load('normali_mano.mat');

% Parametri di analisi
nperm = 2000; % Numero di permutazioni
preclust_pval = 0.05; % P-value per formare cluster (abbassato per maggiore sensibilità)
clust_pval = 0.05; % P-value per la significatività del cluster

ntime_info = 180;  % Limita l'analisi ai primi 180 punti temporali

% Limitare `time_info` ai primi 180 punti
time_info = time_info(1:ntime_info);

% Calcolo della II differenziale (h0) per i dati reali
h0 = zeros(ntime_info, ntime_info);
for t1 = 1:ntime_info
    for t2 = (t1+1):ntime_info
        % Estrazione dei valori di II per tutti i soggetti a t1, t2 nel gruppo ct_ginocchio
        data_ct = cellfun(@(x) x(t1, t2), ct_mano_II);
        % Estrazione dei valori di II per tutti i soggetti a t1, t2 nel gruppo normali_ginocchio
        data_normali = cellfun(@(x) x(t1, t2), normali_mano_II);
        % Esecuzione del t-test tra i due gruppi di dati per la coppia di timepoints
        [~, p, ~, stats] = ttest2(data_ct, data_normali);
        h0(t1, t2) = stats.tstat;
    end
end
h0 = h0 + h0';  % Completa la matrice triangolare superiore

% Calcolo delle permutazioni (hp)
hp = zeros(ntime_info, ntime_info, nperm);
all_II = cat(3, ct_mano_II{:}, normali_mano_II{:});
nsub_ct = numel(ct_mano_II);
nsub_normali = numel(normali_mano_II);
nsub_total = nsub_ct + nsub_normali;

for iperm = 1:nperm
    perm_indices = randperm(nsub_total);
    perm_ct_indices = perm_indices(1:nsub_ct);
    perm_normali_indices = perm_indices(nsub_ct+1:end);
    
    for t1 = 1:ntime_info
        for t2 = (t1+1):ntime_info
            data_perm_ct = all_II(t1, t2, perm_ct_indices);
            data_perm_normali = all_II(t1, t2, perm_normali_indices);
            [~, ~, ~, stats] = ttest2(data_perm_ct, data_perm_normali);
            hp(t1, t2, iperm) = stats.tstat;
        end
    end
    
    % Visualizza la percentuale di completamento
    fprintf('Permutazione %d di %d completata (%.2f%%)\n', iperm, nperm, (iperm/nperm)*100);
end

% Completa la matrice delle permutazioni
for iperm = 1:nperm 
    hp(:,:,iperm) = hp(:,:,iperm) + hp(:,:,iperm)'; 
end 

clust_max = zeros(1, nperm);
for iperm = 1:nperm
    zvals = squeeze((hp(:,:,iperm) - mean(hp(:,:,setdiff(1:nperm, iperm)), 3)) ./ std(hp(:,:,setdiff(1:nperm, iperm)), [], 3));
    zvals(abs(zvals) < norminv(1 - preclust_pval)) = 0;

    clust_info = bwconncomp(zvals);
    if isempty(clust_info.PixelIdxList)
        clust_max(iperm) = 0;
    else
        % Calcola la somma dei valori t all'interno di ogni cluster
        clust_sum = cellfun(@(x) sum(zvals(x)), clust_info.PixelIdxList);
        clust_max(iperm) = max(clust_sum);
    end
end

zmap = (h0 - mean(hp, 3)) ./ std(hp, [], 3);
zmap(abs(zmap) < norminv(1 - preclust_pval)) = 0;

clust_info = bwconncomp(zmap);
clust_size = cellfun(@numel, clust_info.PixelIdxList);
clust_sum = cellfun(@(x) sum(zmap(x)), clust_info.PixelIdxList);

clust_th = prctile(clust_max, 100 - clust_pval * 100);
clust_rem = find(clust_sum < clust_th);
for i = 1:length(clust_rem)
    zmap(clust_info.PixelIdxList{clust_rem(i)}) = 0;
end
zmap(isnan(zmap)) = 0;
zmap = logical(zmap);

%% Plot
figure;
time_plot = time_info / 1000; % Converti il tempo in secondi
imagesc(time_plot, time_plot, h0);
hold on
contour(time_plot, time_plot, zmap, 1, 'linecolor', 'k', 'LineWidth', 1);
colormap(brewermap([], '*RdBu'));
colorbar
axis square
lim = max(abs(min(h0(:))), max(h0(:)));
clim([-lim, lim]) % center zero
xlabel('Time (s)'), ylabel('Time (s)') % Modifica le etichette dell'asse per riflettere il tempo in secondi
title('Cluster-Based Permutation Test Results with p < 0.05');

%% Plot con maschera per colorare solo i cluster significativi
figure;
time_plot = time_info / 1000; % Converti il tempo in secondi

% Applica la maschera: mantieni solo i valori significativi di h0
h0_masked = h0; % Copia di h0
h0_masked(~zmap) = 0; % Imposta a NaN tutto ciò che non è significativo

% Visualizza solo i valori all'interno dei cluster
imagesc(time_plot, time_plot, h0_masked); 
hold on;

% Contorno dei cluster significativi
contour(time_plot, time_plot, zmap, 1, 'linecolor', 'k', 'LineWidth', 1);

% Impostazioni grafiche
colormap(brewermap([], '*RdBu')); % Mantiene la scala di colore
colorbar; % Aggiunge la barra della scala
axis square; % Assicura un aspetto quadrato
caxis([-lim, lim]); % Mantiene il range simmetrico dei valori
xlabel('Time (s)');
ylabel('Time (s)');
title('Cluster-Based Permutation Test Results');

% Imposta lo sfondo bianco
set(gca, 'Color', 'w'); % Imposta lo sfondo bianco per le zone fuori dai cluster

%% clusters info

% Ricalcola la mappa dei cluster dopo la soglia
clust_info_signif = bwconncomp(zmap);

% Calcola la somma dei valori z per i cluster significativi
clust_sum_signif = cellfun(@(x) sum(zmap(x)), clust_info_signif.PixelIdxList);

% Calcola la dimensione di ciascun cluster significativo
clust_size_signif = cellfun(@numel, clust_info_signif.PixelIdxList);

% Numero di cluster significativi
num_signif_clusters = clust_info_signif.NumObjects;

% Prealloca un array di strutture per i cluster significativi
clusters(num_signif_clusters) = struct();

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
    
    % Dimensione del cluster
    cluster_size = numel(cluster_indices);
    
    % Memorizza le informazioni nella struttura
    clusters(i).number = i;
    clusters(i).weight = cluster_weight;
    clusters(i).size = cluster_size;
    clusters(i).t1_min = min(t1_times);
    clusters(i).t1_max = max(t1_times);
    clusters(i).t2_min = min(t2_times);
    clusters(i).t2_max = max(t2_times);
    
    % (Opzionale) Memorizza i tempi completi
    clusters(i).t1_times = t1_times;
    clusters(i).t2_times = t2_times;
    
    % Visualizza le informazioni
    fprintf('Cluster %d:\n', clusters(i).number);
    fprintf(' - Peso (somma dei valori t): %.2f\n', clusters(i).weight);
    fprintf(' - Dimensione: %d punti\n', clusters(i).size);
    fprintf(' - Intervalli temporali coinvolti:\n');
    fprintf('   * Time 1: %.2f ms to %.2f ms\n', clusters(i).t1_min, clusters(i).t1_max);
    fprintf('   * Time 2: %.2f ms to %.2f ms\n\n', clusters(i).t2_min, clusters(i).t2_max);
end

%% Save variables
% Lista delle variabili da salvare
variables_to_save = {'h0', 'zmap', 'time_info', 'lim', 'clust_info_signif', 'clust_sum_signif', 'clust_size_signif', 'clusters'};

% Salva le variabili in un file .mat
save('2kperm_cluster_ct_vs_normali_mano_results.mat', variables_to_save{:});
