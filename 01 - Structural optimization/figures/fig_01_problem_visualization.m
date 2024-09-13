figure('Name', 'Optimization problem visualization');
tiledlayout(1, 2);

% Constraints
nexttile

mesh(X, Y, g1(X, Y), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none');
hold on
mesh(X, Y, g2(X, Y), 'FaceColor', 'g', 'FaceAlpha', 0.5, 'LineStyle', 'none');
mesh(X, Y, g3(X, Y), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none');

axis equal
axis([x([1, end]) y([1, end]) -0.5 0])
title('Constraints curves')
xlabel('x_1')
ylabel('x_2')
legend('$g_1$', '$g_2$', '$g_3$', 'interpreter', 'latex')

% Contours
nexttile
hold on
grid on

imagesc(x, y, f(X, Y));
contour(X, Y, g1(X, Y), [0,0], 'r', 'LineWidth', 2);
contour(X, Y, g2(X, Y), [0,0], 'g', 'LineWidth', 2);
contour(X, Y, g3(X, Y), [0,0], 'b', 'LineWidth', 2);
plot(solution_point.x1, solution_point.x2, 'ok');
text(solution_point.x1 * 1.02, solution_point.x2 * 1.02, sprintf('f(x_1, x_2) = %.3f', solution_point.f))

colorbar

axis('equal');
title('Objective function, constraints and solution');
xlabel('x_1');
ylabel('x_2');
legend('$g_1$', '$g_2$', '$g_3$', 'Solution', 'interpreter', 'latex')