classdef TrussElement
    properties
        node1
        node2
        E % Elastic modulus [N/mm^2]
        A % Cross-sectional area [mm^2]
        L % Length of the element [mm]
        gamma % Element inclination [rad]
    end

    methods
        function obj = TrussElement(node1, node2, E, A)
            % Constructor for TrussElement class
            obj.node1 = node1;
            obj.node2 = node2;
            obj.E = E;
            obj.A = A;
            obj.L = obj.getLength();
            obj.gamma = obj.getGamma();
        end

        function L = getLength(obj)
            % Compute the length of the truss element
            x1 = obj.node1.x; y1 = obj.node1.y;
            x2 = obj.node2.x; y2 = obj.node2.y;
            L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        end

        function gamma = getGamma(obj)
            % Compute the angle (gamma) of the truss element
            x1 = obj.node1.x; y1 = obj.node1.y;
            x2 = obj.node2.x; y2 = obj.node2.y;
            gamma = atan2(y2 - y1, x2 - x1);
        end

        function stiffness = getStiffnessMatrix(obj)
            % Compute the local stiffness matrix for the truss element
            c = cos(obj.gamma);
            s = sin(obj.gamma);
            k = (obj.E * obj.A) / obj.L;

            stiffness = k * [
                 c^2  c*s -c^2 -c*s;
                 c*s  s^2 -c*s -s^2;
                -c^2 -c*s  c^2  c*s;
                -c*s -s^2  c*s  s^2
                ];
        end
    end
end
