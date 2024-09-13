reset(0)
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
% set(0, 'defaultaxesfontsize', 15);
% set(0, 'DefaultLineLineWidth', 2);

plot_struct.flags = true * [1];
% plot_struct.export_path = 'latex/img/MATLAB';
plot_struct.data = cell(0);

if (plot_struct.flags(1))
    run("fig_01_undeformed_structure.m");
end

if (isfield(plot_struct, 'export_path'))
    for plot_idx = 1:numel(plot_struct.data)

        current_plot = plot_struct.data{plot_idx};
        tile = current_plot{1};
        local_path = current_plot{2};

        filename = [plot_struct.export_path local_path '.png'];
        exportgraphics(tile, filename, 'Resolution', 300);

    end
end

clear plot_struct plot_idx current_plot local_path filename tile 