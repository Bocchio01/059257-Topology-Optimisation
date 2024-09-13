clc
clear variables
close all

%% Structure and Applied loads

[structure_brk, structure_acc] = deal(load_structure('models/motorcycle_frame.inp'));

% Braking
structure_brk.nodes(8).forceX = -13000;
structure_brk.nodes(8).forceY = +2000;
structure_brk.nodes(9).forceX = +9000;
structure_brk.nodes(9).forceY = +0;

% Acceleration
structure_acc.nodes(1).forceX = +2000;
structure_acc.nodes(1).forceY = -8000;


%% Optimization problem constraints

alpha = 0.1; % []
lower_bound = 1e-3; % [mm^2]
upper_bound = 700;  % [mm^2]
w1 = 0.5;
w2 = 1 - w1;

max_volume = upper_bound * sum([structure_brk.elements.L], 'all');
max_volume_constrained = max_volume * alpha;


%% Problem solution

epsilon = 1e-5;
N = 10;
lambda = 1;
results = struct( ...
    'C', zeros(N, 1), ...
    'V', zeros(N, 1), ...
    'alpha', zeros(N, 1), ...
    'dCdA', zeros(N, length(structure_brk.elements)), ...
    'A', zeros(N, length(structure_brk.elements)), ...
    'stress', zeros(N, length(structure_brk.elements)));

for k = 1 : N

    % FEM Analysis
    structure_brk = structure_brk.runFEM();
    structure_acc = structure_acc.runFEM();
    results.C(k) = w1 * structure_brk.F' * structure_brk.U + ...
                   w2 * structure_acc.F' * structure_acc.U;

    % Sensitivity Analysis
    for ii = 1:length(structure_brk.elements)

        k01 = structure_brk.elements(ii).getStiffnessMatrix / structure_brk.elements(ii).A;
        ue1 = structure_brk.U(structure_brk.getDOFidxs(structure_brk.elements(ii)));

        k02 = structure_acc.elements(ii).getStiffnessMatrix / structure_acc.elements(ii).A;
        ue2 = structure_acc.U(structure_acc.getDOFidxs(structure_acc.elements(ii)));

        results.dCdA(k, ii) = w1 * (- ue1' * k01 * ue1) + ...
                              w2 * (- ue2' * k02 * ue2);

    end

    % Dual problem solution
    handler = @(lambda) dphi( ...
        lambda, ...
        [structure_brk.elements.A]', ...
        lower_bound, ...
        upper_bound, ...
        results.dCdA(k, :)', ...
        [structure_brk.elements.L]', ...
        max_volume_constrained);
    
    % Solve sum(L * A(lambda*)) - V_{0} == 0 and compute A(lambda*)
    lambda = fsolve(handler, lambda);
    [~, A] = handler(lambda);

    % Variables updates
    for ii = 1:length(structure_brk.elements)
        structure_brk.elements(ii).A = A(ii);
        structure_acc.elements(ii).A = A(ii);
    end
    results.V(k) = A' * [structure_brk.elements.L]';
    results.alpha(k) = results.V(k) / max_volume;
    results.A(k, :) = A;
    % results.stress(k, :) = structure1.computeStresses();

    % Check convergence condition
    if(norm(results.dCdA(k, :)) < epsilon)
        break;
    end

end

clear ii k
clear epsilon N A handler lambda
clear ue k0 dofIndices node2Index node1Index


%% Plots

structure = structure_brk;
run("figures\routine.m")


%% Functions

function [y, x] = dphi(lambda, A, lower_bound, upper_bound, dCdA, l_elements, max_volume_constrained)

lambda(lambda < 0) = 0;

q0 = -A.^2 .* dCdA;
x = sqrt(q0 ./ (lambda * l_elements));

x(x > upper_bound) = upper_bound;
x(x < lower_bound) = lower_bound;

y = l_elements' * x - max_volume_constrained;

end
