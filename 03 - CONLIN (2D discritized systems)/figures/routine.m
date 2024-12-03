reset(0)
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 15);
set(0, 'DefaultLineLineWidth', 2);

% plot_struct.export_path = 'latex/img/MATLAB';
plot_struct.data = cell(0);

% Undeformed
% figure_undeformed_structure = figure('Name', 'Undeformed structure');
% 
% nexttile
% hold on
% grid on
% 
% plot_structure(undeformed_structure.runFEM(), struct('color', '#D95319', 'LineWidth', 2, 'force', true));
% 
% axis equal
% axis padded
% title('Undeformed structure')
% xlabel('[mm]')
% ylabel('[mm]')
% ylim([-500 4500])
% xlim([-500 8500])


% Results
figure_topology_results = figure('Name', 'Topology optimization results');
tile = tiledlayout(2, 3);

% Compliance
nexttile(tile, 1)
hold on
grid on

plot(results.C(1:k), 'ob');

title(['Optimization compliance @alpha=' num2str(alpha, '%.2f')])
xlabel('Iteration #')
ylabel('Compliance [Nmm]')

% Volume
nexttile(tile, 4)
hold on
grid on

plot(results.alpha(1:k) * 100, 'ob');

ylim([0 30])

title(['Volume fraction @alpha=' num2str(alpha, '%.2f')])
xlabel('Iteration #')
ylabel('Volume fraction []')

ytickformat('percentage')


% Deformed
nexttile(tile, 2, [2, 2])
hold on
grid on

plot_structure(structure);

title(['Deformed structure @alpha=' num2str(alpha, '%.2f')])
xlabel('[mm]')
ylabel('[mm]')


% plot_struct.data{end+1} = {figure_undeformed_structure, '/undeformed_structure'};
plot_struct.data{end+1} = {figure_topology_results, ['/results_alpha' num2str(alpha*100)]};


if (isfield(plot_struct, 'export_path'))
    for plot_idx = 1:numel(plot_struct.data)

        current_plot = plot_struct.data{plot_idx};
        tile = current_plot{1};
        local_path = current_plot{2};

        filename = [plot_struct.export_path local_path '.pdf'];
        exportgraphics(tile, filename, 'ContentType', 'vector');

    end
end

clear plot_struct plot_idx current_plot local_path filename tile 