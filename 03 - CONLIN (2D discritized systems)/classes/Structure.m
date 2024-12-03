classdef Structure
% Structure class for Finite Element Analysis (FEA) of truss systems
%
% This class represents a 2D truss structure, consisting of nodes and
% elements (truss members). It is capable of assembling the global
% stiffness matrix, applying constraints, and computing displacements
% and stresses for each element using the Finite Element Method (FEM).
%
% Example:
%   nodes       = [Node(...), Node(...)];
%   elements    = [Element(...), Element(...)];
%   structure   = Structure(nodes, elements);
%   structure.runFEM();
%   stresses    = structure.computeStresses();


    properties
        nodes % Array of Node objects representing the structure's joints.
        elements % Array of Element objects representing truss members.
        K % Global stiffness matrix (2 DOF per node).
        F % Global load vector (forces applied at nodes).
        U % Global displacement vector (displacements at nodes).
    end

    methods
        function obj = Structure(nodes, elements)
            % Constructor, initializes the nodes and elements of the
            % structure

            obj.nodes = nodes;
            obj.elements = elements;
        end

        function obj = runFEM(obj)
            % Runs the full finite element analysis, including assembling
            % stiffness matrix, load vector, applying constraints, and
            % computing displacements

            obj.K = obj.assembleGlobalStiffness();
            obj.F = obj.assembleLoadVector();
            obj.applyConstraints();
            obj.U = obj.computeDisplacements();
        end

        function K = assembleGlobalStiffness(obj)
            % Assembles the global stiffness matrix based on the
            % individual element stiffness matrices

            K = zeros(2 * length(obj.nodes));

            for i = 1:length(obj.elements)

                dofIndices = obj.getDOFidxs(obj.elements(i));
                K(dofIndices, dofIndices) = K(dofIndices, dofIndices) + obj.elements(i).getStiffnessMatrix();
            
            end

        end

        function F = assembleLoadVector(obj)
            % Assembles the global load vector from the forces 
            % applied to each node.
            
            F = zeros(2 * length(obj.nodes), 1);
            F(1:2:end) = [obj.nodes.forceX];
            F(2:2:end) = [obj.nodes.forceY];

        end

        function applyConstraints(obj)
            % Applies boundary conditions (constraints) to the global 
            % stiffness matrix and load vector

            dofIndex = [
                find([obj.nodes.isConstrainedX] == true) * 2 - 1, ...
                find([obj.nodes.isConstrainedY] == true) * 2
            ];

            obj.K(dofIndex, :) = 0;
            obj.K(:, dofIndex) = 0;
            obj.K(dofIndex, dofIndex) = eye(length(dofIndex));
            obj.F(dofIndex) = 0;

        end

        function U = computeDisplacements(obj)
            % Solves for the global displacement vector using the 
            % global stiffness matrix and load vector

            U = zeros(size(obj.K, 1), 1);
            freeDOFs = ~logical([ ...
                obj.nodes.isConstrainedX; ...
                obj.nodes.isConstrainedY ...
                ]);
            K_FF = obj.K(freeDOFs, freeDOFs);
            F_F = obj.F(freeDOFs);

            U(freeDOFs) = K_FF \ F_F;
        end

        function sigma = computeStresses(obj)
            % Computes the stress in each element of the truss 
            % based on the displacement results

            sigma = zeros(length(obj.elements), 1);

            for ii = 1:length(obj.elements)
                elem = obj.elements(ii);

                ev = [cos(elem.gamma) sin(elem.gamma)];
                B = [-ev ev];

                DOFidxs = obj.getDOFidxs(elem);
                uel = obj.U(DOFidxs);

                sigma(ii) = (elem.E / elem.L) * B * uel;
            end
        end

        function DOFidxs = getDOFidxs(obj, elements)

            DOFidxs = NaN;
   
            for ii = 1 : length(elements)

            node1Index = find(obj.nodes == elements(ii).node1);
            node2Index = find(obj.nodes == elements(ii).node2);

            DOFidxs = [DOFidxs, 2 * node1Index - 1, 2 * node1Index, 2 * node2Index - 1, 2 * node2Index];

            end

            DOFidxs = unique(DOFidxs(~isnan(DOFidxs)), 'stable');
        end
    end
end
