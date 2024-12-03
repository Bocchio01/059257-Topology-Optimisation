function plot_structure(structure, opts)

if (nargin == 1)
    opts = struct('color', false, 'LineWidth', false, 'force', true);
end

% Beams
max_A = max([structure.elements.A]);
sigma = structure.computeStresses();
for ii = 1:length(structure.elements)

    element = structure.elements(ii);

    if sigma(ii) > 0
        color = 'r';
    else
        color = 'b';
    end

    if (opts.color)
        color = opts.color;
    end

    LineWidth = element.A / 50;
    if (opts.LineWidth)
        LineWidth = opts.LineWidth;
    end

    plot([element.node1.x element.node2.x], ...
        [element.node1.y element.node2.y], ...
        'LineWidth', LineWidth, ...
        'Color', color);

end

% Graphical constraints
radius = 3;
dimension = 15;
shape_horizontal_constraint = [ 0 0; -sqrt(3)/2 .5; -sqrt(3)/2 -.5; 0 0] * dimension;
shape_vertical_constraint = [ 0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0] * dimension;
triangolo_r = [0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0] * dimension * [sqrt(2)/2 -sqrt(2/2); -sqrt(2)/2 -sqrt(2)/2];

% Nodes
for ii = 1:length(structure.nodes)

    node = structure.nodes(ii);

    fplot(@(t) radius*sin(t) + node.x, @(t) radius*cos(t) + node.y, 'Color', 'k', 'LineWidth', 2);
    text(node.x + 5*radius, node.y + 5*radius, num2str(ii), 'Color', 'k');

    if(logical(node.isConstrainedX))
        fill(node.x + shape_horizontal_constraint(:, 1), node.y + shape_horizontal_constraint(:, 2), 'k');
    end

    if(logical(node.isConstrainedY))
        fill(node.x + shape_vertical_constraint(:, 1), node.y + shape_vertical_constraint(:, 2), 'k');
    end

    if(opts.force && logical(node.forceX))
        q = quiver(node.x, node.y, 3*dimension * sign(node.forceX), 0, 'k', 'LineWidth', 2);
        q.MaxHeadSize = 1;
    end

    if(opts.force && logical(node.forceY))
        q = quiver(node.x, node.y, 0, 3*dimension * sign(node.forceY), 'k', 'LineWidth', 2);
        q.MaxHeadSize = 1;
    end

end

axis padded

end