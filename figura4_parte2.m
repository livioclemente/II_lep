clear; clc;

% Carica i dati dei gruppi
load('ct_ginocchio.mat');
load('dl_pth_ginocchio.mat');

% Imposta il canale e il timepoint da visualizzare
chi = 11; % canale specifico
time_ms = 350; % tempo in millisecondi
[~, ti] = min(abs(time_info - time_ms)); % trova l'indice del timepoint più vicino a 400 ms

time_window = [350, 450]; % finestra temporale in millisecondi
window_indices = find(time_info >= time_window(1) & time_info <= time_window(2));

% Trova il numero minimo di soggetti tra i due gruppi
min_subjects = min(length(ct_ginocchio_LEP), length(dl_pth_ginocchio_LEP));

% Numero di soggetti da visualizzare (non superiore a min_subjects)
Nplt = min(10, min_subjects);
cols = {'r' 'b'};

% Visualizzazione dei tracciati EEG
figure;
hold on;
offset = 4;

% Seleziona soggetti casuali da ciascun gruppo
ct_idx = randperm(length(ct_ginocchio_LEP), Nplt);
dl_pth_idx = randperm(length(dl_pth_ginocchio_LEP), Nplt);

% Visualizza i tracciati EEG per il gruppo 'ct'
for pii = 1:Nplt
    pi = ct_idx(pii);
    eeg_ct = ct_ginocchio_LEP{pi};
    plot(time_info / 1000, (pii - 1) * offset + eeg_ct(chi, :), 'r'); % colore rosso per ct
    plot(time_info(ti) / 1000, (pii - 1) * offset + eeg_ct(chi, ti), '.', 'MarkerSize', 16, 'Color', 'r');
end

% Visualizza i tracciati EEG per il gruppo 'dl_pth'
for pii = 1:Nplt
    pi = dl_pth_idx(pii);
    eeg_dl_pth = dl_pth_ginocchio_LEP{pi};
    plot(time_info / 1000, (pii - 1) * offset + eeg_dl_pth(chi, :), 'b'); % colore blu per dl_pth
    plot(time_info(ti) / 1000, (pii - 1) * offset + eeg_dl_pth(chi, ti), '.', 'MarkerSize', 16, 'Color', 'b');
end

xlim([time_info(1) / 1000, time_info(end) / 1000]); % in secondi
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.025]);

% Aggiungi una legenda corretta
h_ct = plot(nan, nan,  'bo', 'MarkerSize', 8, 'MarkerFaceColor','r'); % Placeholder per linea rossa
h_dl_pth = plot(nan, nan, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Placeholder per pallino blu
legend([h_ct, h_dl_pth], {'CT', 'DLPTH'}, 'Location', 'best');

% Aggiungi un riquadro attorno ai punti di interesse a 400 ms
rectangle('Position', [time_info(ti)/1000 - 0.025, -offset, 0.05, offset * 2 * Nplt], 'EdgeColor', 'k', 'LineWidth', 1.5);

title(sprintf('EEG Traces at %d ms', time_ms));
xlabel('Time (s)');
ylabel('EEG Amplitude');

%%

% Estrazione degli ensemble per il canale specifico e la finestra temporale specifica
ens_ct = [];
for pii = 1:Nplt
    pi = ct_idx(pii);
    eeg_ct = ct_ginocchio_LEP{pi};
    ens_ct = [ens_ct, eeg_ct(chi, ti)];
end

ens_dl_pth = [];
for pii = 1:Nplt
    pi = dl_pth_idx(pii);
    eeg_dl_pth = dl_pth_ginocchio_LEP{pi};
    ens_dl_pth = [ens_dl_pth, eeg_dl_pth(chi, ti)];
end

% Calcola la densità di probabilità
kxi = linspace(-3, 3, 100);
[k_ct, ~] = ksdensity(ens_ct(:), kxi);
[k_dl_pth, ~] = ksdensity(ens_dl_pth(:), kxi);

% Visualizzazione della densità di probabilità
figure;
plot(kxi, k_ct, 'r', 'LineWidth', 2);
hold on;
plot(kxi, k_dl_pth, 'b', 'LineWidth', 2);
box off;
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.025]);
title('Densità di Probabilità');
xlabel('Valori EEG');
ylabel('Densità di Probabilità');
legend({'CT', 'DLPTH'}, 'Location', 'best');

