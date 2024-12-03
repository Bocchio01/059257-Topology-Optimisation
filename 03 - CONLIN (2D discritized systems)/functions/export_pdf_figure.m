function export_pdf_figure(plot_struct)

if (isfield(plot_struct, 'export_path'))
    for plot_idx = 1:numel(plot_struct.data)

        current_plot = plot_struct.data{plot_idx};
        tile = current_plot{1};
        local_path = current_plot{2};

        filename = [plot_struct.export_path local_path '.pdf'];
        exportgraphics(tile, filename, 'ContentType', 'vector');

    end
end

end
