clear; clc;

addpath('/Users/livioclemente/Documents/MATLAB/bluewhitered/');
addpath('/Users/livioclemente/Documents/MATLAB/eeglab2023.0/');
addpath('/Users/livioclemente/Documents/MATLAB/EKG_EEG-master/');
addpath('/Users/livioclemente/Documents/MATLAB/DrosteEffect_BrewerMap_3/');

% Carica i dati dei gruppi
load('ct_ginocchio.mat');
load('normali_ginocchio.mat');

% Parametri di analisi
nperm = 200; % Numero di permutazioni
preclust_pval = 0.2; % P-value per formare cluster (abbassato per maggiore sensibilità)
clust_pval = 0.2; % P-value per la significatività del cluster

% Ottenere le dimensioni per le matrici di II
ntime_info = size(ct_ginocchio_II{1}, 1);  % Assume tutte le matrici II hanno la stessa dimensione

% Calcolo della II differenziale (h0) per i dati reali
h0 = zeros(ntime_info, ntime_info);
for t1 = 1:ntime_info
    for t2 = (t1+1):ntime_info
        % Estrazione dei valori di II per tutti i soggetti a t1, t2 nel gruppo ct_ginocchio
        data_ct = cellfun(@(x) x(t1, t2), ct_ginocchio_II);
        % Estrazione dei valori di II per tutti i soggetti a t1, t2 nel gruppo normali_ginocchio
        data_normali = cellfun(@(x) x(t1, t2), normali_ginocchio_II);
        % Esecuzione del t-test tra i due gruppi di dati per la coppia di timepoints
        [~, p, ~, stats] = ttest2(data_ct, data_normali);
        h0(t1, t2) = stats.tstat;
    end
end
h0 = h0 + h0';  % Completa la matrice triangolare superiore

% Calcolo delle permutazioni (hp)
hp = zeros(ntime_info, ntime_info, nperm);
all_II = cat(3, ct_ginocchio_II{:}, normali_ginocchio_II{:});
nsub_ct = numel(ct_ginocchio_II);
nsub_normali = numel(normali_ginocchio_II);
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

% Plot
time_plot = time_info / 1000; % Converti il tempo in secondi
imagesc(time_plot, time_plot, h0);
hold on
contour(time_plot, time_plot, zmap, 1, 'linecolor', 'k', 'LineWidth', 1);
colormap(brewermap([], '*RdBu'));
colorbar
axis square
lim = max(abs(min(h0(:))), max(h0(:)));
caxis([-lim, lim]) % center zero
xlabel('Time (s)'), ylabel('Time (s)') % Modifica le etichette dell'asse per riflettere il tempo in secondi
title('Cluster-Based Permutation Test Results with p < 0.05');

