clear; clc; close all;

addpath('/Users/livioclemente/Documents/MATLAB/bluewhitered/');
addpath('/Users/livioclemente/Documents/MATLAB/eeglab2023.0/');
addpath('/Users/livioclemente/Documents/MATLAB/EKG_EEG-master/');
addpath('/Users/livioclemente/Documents/MATLAB/DrosteEffect_BrewerMap_3/');

time_start = -0.1; % start time in seconds
time_end = 1; % end time in seconds
srate = 256;
time_epoch = time_start:1/srate:time_end-1/srate;
n_samples = round((time_end - time_start) * srate);
time_info = linspace(time_start * 1000, time_end * 1000 - 1 / srate * 1000, n_samples);
ntime_info = length(time_info);

% Tracciati per il ginocchio
data_path_knee_ct = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/MI_LEP/conditions/ct/ginocchio/';
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

for isub = 1:nsub_knee_ct % Carica i dati
    EEG_knee_ct = pop_loadset(sub_knee_ct(isub).name, data_path_knee_ct);
    DATA_knee_ct{isub} = EEG_knee_ct.data;

    [electrodes, Nt, trials_knee] = size(DATA_knee_ct{isub}); % Estrazione e normalizzazione dei dati
    x = squeeze(DATA_knee_ct{isub}(11, :, :))';
    vas = repmat(vas_values(isub), trials_knee, 1); % Ripeti il valore VAS per ogni trial per il soggetto corrente
    x_knee_ct = copnorm(x);
    vas_knee_ct = copnorm(vas);

    x_values_knee_ct{isub} = x_knee_ct;
    vas_values_knee_ct{isub} = vas_knee_ct;

    MI = zeros(1, Nt);
    for ti = 1:Nt
        MI(ti) = mi_gg(x_knee_ct(:, ti), vas_knee_ct(:, 1), true, true);
    end
    
    MI_knee_ct{isub} = MI;

    noise = .00000005 * randn(size(x_knee_ct, 1), 1);
    II = zeros(Nt, Nt);
    for t1 = 1:Nt
        for t2 = (t1 + 1):Nt
            JMI = mi_gg([x_knee_ct(:, t1) x_knee_ct(:, t2) + noise], vas_knee_ct(:, 1), true, true);
            II(t1, t2) = JMI - MI(t1) - MI(t2);
        end
    end
    II = II + II';
    II_knee_ct{isub} = II;
    
    % LEP
    LEP_knee_ct{isub} = mean(DATA_knee_ct{isub}, 3);
end

%% Global

MI_avg_knee_ct = mean(cat(1, MI_knee_ct{:}), 1);
II_avg_knee_ct = mean(cat(3, II_knee_ct{:}), 3);
LEP_avg_knee_ct = mean(cat(3, LEP_knee_ct{:}), 3);

%% Calcolo Aggregato e Permutazioni

% Preallocazione per aggregazione
all_x_knee_ct = [];
all_vas_knee_ct = [];

% Caricamento e aggregazione dei dati
for isub = 1:nsub_knee_ct
    EEG_knee_ct = pop_loadset(sub_knee_ct(isub).name, data_path_knee_ct);
    x = squeeze(EEG_knee_ct.data(11, :, :))'; % Estrazione dati per canale specifico
    vas = repmat(vas_values(isub), size(x, 1), 1); % Ripeti il valore VAS per ogni trial
    
    % Normalizzazione e aggregazione
    all_x_knee_ct = [all_x_knee_ct; copnorm(x)]; % Aggiungi i dati normalizzati
    all_vas_knee_ct = [all_vas_knee_ct; copnorm(vas)]; % Aggiungi stimoli normalizzati
end

% Calcolo della MI aggregata
MI_aggregated = zeros(1, size(all_x_knee_ct, 2));
for ti = 1:size(all_x_knee_ct, 2)
    MI_aggregated(ti) = mi_gg(all_x_knee_ct(:, ti), all_vas_knee_ct, true, true);
end

nperm = 300; % number of permutations
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
time_plot = time_info*1000;
imagesc(time_plot,time_plot,h0(:,:))
hold on
contour(time_plot,time_plot,zmap,1,'linecolor','k','LineWidth',1)
colormap(brewermap([],'*RdBu'));
colorbar
axis square
lim = max(abs(min(min(h0(:,:)))),max(max(h0(:,:))));
clim([-lim,lim]) % center zero
xlabel('Time (ms)'),ylabel('Time (ms)')
