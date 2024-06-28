%% Plot della MI

load figura_4_humanbrainmapping/gruppi_fig4/ct_ginocchio.mat
load figura_4_humanbrainmapping/gruppi_fig4/dl_pth_ginocchio.mat

% Visualizzazione
figure;

% Tracciare la mutua informazione media per il gruppo 'ct'
subplot(3, 1, 1);
imagesc(time_info, 1:64, ct_ginocchio_MI_all_channels_avg');
title('Mutua Informazione Media - Gruppo CT');
xlabel('Tempo (ms)');
ylabel('Canali');
colorbar;
set(gca, 'XLim', [time_info(1), time_info(end)]); % Imposta i limiti dell'asse x
set(gca, 'YLim', [1 64]); % Imposta i limiti dell'asse y

% Tracciare la mutua informazione media per il gruppo 'dl_pth'
subplot(3, 1, 2);
imagesc(time_info, 1:64, dl_pth_ginocchio_MI_all_channels_avg');
title('Mutua Informazione Media - Gruppo DLPTH');
xlabel('Tempo (ms)');
ylabel('Canali');
colorbar;
set(gca, 'XLim', [time_info(1), time_info(end)]); % Imposta i limiti dell'asse x
set(gca, 'YLim', [1 64]); % Imposta i limiti dell'asse y

% Tracciare la differenza di mutua informazione media tra i due gruppi
subplot(3, 1, 3);
imagesc(time_info, 1:64, (ct_ginocchio_MI_all_channels_avg - dl_pth_ginocchio_MI_all_channels_avg)');
title('Differenza di Mutua Informazione Media (CT - DLPTH)');
xlabel('Tempo (ms)');
ylabel('Canali');
colorbar;
set(gca, 'XLim', [time_info(1), time_info(end)]); % Imposta i limiti dell'asse x
set(gca, 'YLim', [1 64]); % Imposta i limiti dell'asse y

%% Video della MI

% Impostazione dei limiti per la mappa
minlim = min(min(ct_ginocchio_MI_all_channels_avg(:)), min(dl_pth_ginocchio_MI_all_channels_avg(:)));
maxlim = max(max(ct_ginocchio_MI_all_channels_avg(:)), max(dl_pth_ginocchio_MI_all_channels_avg(:)));

% Creazione del video per il gruppo 'ct'
figure;
clear F1;
for ti = 1:length(time_info)
    cla;
    topoplot(ct_ginocchio_MI_all_channels_avg(ti, :), chanlocs, 'maplimits', [minlim maxlim]);
    title(sprintf('%d ms - Gruppo CT', time_info(ti)), 'FontSize', 12);
    colorbar;
    F1(ti) = getframe(gcf);
end

fname1 = 'ct_ginocchio_info.mov';
v1 = VideoWriter(fname1, 'MPEG-4');
v1.FrameRate = 5;
open(v1);
for i = 1:length(F1)
    writeVideo(v1, F1(i));
end
close(v1);
close(gcf);

% Creazione del video per il gruppo 'dl_pth'
figure;
clear F2;
for ti = 1:length(time_info)
    cla;
    topoplot(dl_pth_ginocchio_MI_all_channels_avg(ti, :), chanlocs, 'maplimits', [minlim maxlim]);
    title(sprintf('%d ms - Gruppo DL_PTH', time_info(ti)), 'FontSize', 12);
    colorbar;
    F2(ti) = getframe(gcf);
end

fname2 = 'dl_pth_ginocchio_info.mov';
v2 = VideoWriter(fname2, 'MPEG-4');
v2.FrameRate = 5;
open(v2);
for i = 1:length(F2)
    writeVideo(v2, F2(i));
end
close(v2);
close(gcf);

%% Plot ERP

% Impostazione dei limiti per la mappa
minlim = min(min(ct_ginocchio_MI_all_channels_avg(:)), min(dl_pth_ginocchio_MI_all_channels_avg(:)));
maxlim = max(max(ct_ginocchio_MI_all_channels_avg(:)), max(dl_pth_ginocchio_MI_all_channels_avg(:)));

% Visualizzazione degli ERP Mediati per un Canale Specifico
ci = 11; % canale specifico
time = time_info / 1000; % converti time_info in secondi
tidx = find((time > -0.1) & (time < 1)); % periodo che vuoi visualizzare

% Media degli ERP per ciascun gruppo
ct_erp_mean = ct_ginocchio_LEP_avg(ci, :);
dl_pth_erp_mean = dl_pth_ginocchio_LEP_avg(ci, :);

% Visualizzazione degli ERP Mediati
figure;
subplot(2, 1, 1);
hold on;
plot(time(tidx), ct_erp_mean(tidx), 'r');
plot(time(tidx), dl_pth_erp_mean(tidx), 'b');
legend('CT', 'DLPTH');
title('ERP Mediati - Gruppo CT vs DLPTH');

subplot(2, 1, 2);
hold on;
plot(time(tidx), dl_pth_erp_mean(tidx) - ct_erp_mean(tidx), 'k');
legend('Differenza ERP DLPTH - CT');
title('Differenza di ERP');

% Visualizzazione della MI
figure;
subplot(3, 1, 1);
plot(time_info, ct_ginocchio_MI_all_channels_avg(:, ci));
title('MI - Gruppo CT');
xlabel('Tempo (ms)');
ylabel('MI');
set(gca, 'XLim', [time_info(1), time_info(end)]);

subplot(3, 1, 2);
plot(time_info, dl_pth_ginocchio_MI_all_channels_avg(:, ci));
title('MI - Gruppo DL_PTH');
xlabel('Tempo (ms)');
ylabel('MI');
set(gca, 'XLim', [time_info(1), time_info(end)]);

subplot(3, 1, 3);
plot(time_info, ct_ginocchio_MI_all_channels_avg(:, ci) - dl_pth_ginocchio_MI_all_channels_avg(:, ci));
title('Differenza di MI (CT - DL_PTH)');
xlabel('Tempo (ms)');
ylabel('MI');
set(gca, 'XLim', [time_info(1), time_info(end)]);


