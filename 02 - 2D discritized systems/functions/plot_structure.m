function plot_structure(structure)

% Beams
for ii = 1:length(structure.elements)

    element = structure.elements(ii);

    plot([element.node1.x element.node2.x], ...
        [element.node1.y element.node2.y], ...
        'LineWidth', element.A / 1e2, ...
        'Color', 'k');

end

% Nodes
for ii = 1:length(structure.nodes)

    radius = 3;
    node = structure.nodes(ii);

    fplot(@(t) radius*sin(t) + node.x, @(t) radius*cos(t) + node.y, 'Color', 'r', 'LineWidth', 2);
    text(node.x + radius, node.y + radius, num2str(ii), 'Color', 'r');

    if(logical(node.isConstrainedX))
        plot(nsidedpoly(3, 'Center', [node.x, node.y], 'SideLength', 2 * radius));
    end

    if(logical(node.isConstrainedY))
        plot(nsidedpoly(3, 'Center', [node.x, node.y], 'SideLength', 2 * radius));
    end

end

axis equal

end