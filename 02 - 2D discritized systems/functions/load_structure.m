function [structure] = load_structure(filename)

lines = readlines(filename);

% NODES
line_idx = find(contains(lines, "*NODES", 'IgnoreCase', true), 1, 'first') + 2;
while (~contains(lines(line_idx), ["*ENDNODES", "!"], 'IgnoreCase', true))

    tmp = sscanf(replace(lines(line_idx), ',', '.'), '%i %i %i %f %f');

    nodes(tmp(1)) = Node(tmp(4), tmp(5), tmp(2), tmp(3));

    line_idx = line_idx + 1;

end

% PROPERTIES
line_idx = find(contains(lines, "*PROPERTIES", 'IgnoreCase', true), 1, 'first') + 2;
while (~contains(lines(line_idx), "*ENDPROPERTIES", 'IgnoreCase', true))

    tmp = sscanf(lines(line_idx), '%i %f %f');
    properties(tmp(1), :) = tmp(2:3);

    line_idx = line_idx + 1;

end

% BEAMS
line_idx = find(contains(lines, "*BEAMS", 'IgnoreCase', true), 1, 'first') + 2;
while (~contains(lines(line_idx), "*ENDBEAMS", 'IgnoreCase', true))

    tmp = sscanf(lines(line_idx), '%i %i %i %i');

    E = properties(tmp(4), 1);
    A = properties(tmp(4), 2);
    elements(tmp(1)) = TrussElement(nodes(tmp(2)), nodes(tmp(3)), E, A);
    
    line_idx = line_idx + 1;

end

structure = Structure(nodes, elements);

end