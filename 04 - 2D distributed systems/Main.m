clc
clear variables
close all

%% Plots setup

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 15);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultFigureColormap', gray)


%% First requests

figure_p = figure('Name', 'Changing p');
tiledlayout(3, 1)

for p = [1.0 2.0 3.0]
    
    rmin = 1.5;
    nexttile
    imagesc(-TO1_distributed(60, 20, 0.5, p, rmin))
    title(['p=' num2str(p) ' | ' 'r_{min}=' num2str(rmin)])
    axis equal
    axis tight

end


figure_rmin = figure('Name', 'Changing rmin');
tiledlayout(3, 1)

for rmin = [1.0 1.2 1.5]
    
    p = 3.0;
    nexttile
    imagesc(-TO1_distributed(60, 20, 0.5, p, rmin))
    title(['p=' num2str(p) ' | ' 'r_{min}=' num2str(rmin)])
    axis equal
    axis tight

end


%% Second request

figure_load_case_2 = figure('Name', 'Fixed left side and extremity force');
tiledlayout(2, 1)

nexttile
imagesc(-TO2_distributed(32, 20, 0.4, 3.0, 1.2))
title(['p=' num2str(3.0) ' | ' 'r_{min}=' num2str(1.2)])
axis equal
axis tight


%% Third request

figure_load_case_3 = figure('Name', 'Fixed left side and extremity force');
tiledlayout(2, 4)

nexttile
imagesc(-TO31_distributed(30, 30, 0.4, 3.0, 1.2))
title(['p=' num2str(3.0) ' | ' 'r_{min}=' num2str(1.2)])
axis equal
axis tight

nexttile
imagesc(-TO32_distributed(30, 30, 0.4, 3.0, 1.2))
title(['p=' num2str(3.0) ' | ' 'r_{min}=' num2str(1.2)])
axis equal
axis tight


%% Fourth request

figure_load_case_4 = figure('Name', 'Fixed left side and extremity force');
tiledlayout(2, 1)

nexttile
imagesc(-TO4_distributed(45, 30, 0.5, 3.0, 1.5))
title(['p=' num2str(3.0) ' | ' 'r_{min}=' num2str(1.5)])
axis equal
axis tight


%% Exports

plot_struct.export_path = 'latex/img/MATLAB';
plot_struct.data = cell(0);

% plot_struct.data{end+1} = {figure_p, '/req_1p'};
% plot_struct.data{end+1} = {figure_rmin, '/req_1rmin'};
% plot_struct.data{end+1} = {figure_load_case_2, '/req_2'};
plot_struct.data{end+1} = {figure_load_case_3, '/req_3'};
% plot_struct.data{end+1} = {figure_load_case_4, '/req_4'};

export_pdf_figure(plot_struct);
