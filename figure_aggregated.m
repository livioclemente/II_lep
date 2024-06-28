clear; clc; close all;

% Definizione dei gruppi e dei distretti
groups = {'ct', 'dl_pth', 'normali', 'pth'};
districts = {'ginocchio', 'mano', 'piede'};

% Impostazioni per i plot
lw = 1.5; % Spessore delle linee nei plot

% Ciclo attraverso ogni gruppo e distretto
for g = 1:length(groups)
    for d = 1:length(districts)
        group = groups{g};
        district = districts{d};
        file_name = sprintf('%s_%s.mat', group, district);
        load(file_name); % Assume the variables inside are named consistently like in the save command

        % Crea una nuova figura
        figure;
        sgtitle(sprintf('Group: %s, District: %s', group, district)); % Titolo generale per la figura

        % Subplot principale: Interaction Information
        axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
        eval(sprintf('imagesc(time_info, time_info, %s_II_aggregated);', sprintf('%s_%s', group, district)));
        clim_max = max(abs(eval(sprintf('%s_II_aggregated(:)', sprintf('%s_%s', group, district)))));
        clim([-clim_max, clim_max]);
        colormap(axm, bluewhitered);
        colorbar;
        title('Interaction Information');

        % Secondo subplot: Average ERP of all conditions
        ax2 = subplot(5, 5, [22 23 24 24.55]);
        eval(sprintf('plot(time_info, %s_LEP_avg(11, :), ''k'', ''LineWidth'', lw);', sprintf('%s_%s', group, district)));
        axis tight;
        xlim([min(time_info), max(time_info)]);
        xlabel('Time (ms)');
        ylabel('LEP');
        box off;

        % Terzo subplot: Mutual Information as a function of time
        ax3 = subplot(5, 5, [1 6 11 16]);
        eval(sprintf('plot(time_info, %s_MI_aggregated, ''k'', ''LineWidth'', lw);', sprintf('%s_%s', group, district)));
        axis tight;
        box off;
        xlabel('Time (ms)');
        ylabel('MI');
        set(gca, 'CameraUpVector', [-1 0 0]);

        % Salva la figura
        save_fig_name = sprintf('%s_%s_figure.png', group, district);
        saveas(gcf, save_fig_name);
        close(gcf); % Chiudi la figura per evitare l'accumulo di figure aperte
    end
end
