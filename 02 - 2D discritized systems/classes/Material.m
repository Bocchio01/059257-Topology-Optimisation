classdef Material
    properties
        E % Elastic modulus (Young's modulus)
        nu % Poisson's ratio
    end
    
    methods
        function obj = Material(E, nu)
            % Constructor for Material class
            obj.E = E;
            obj.nu = nu;
        end
    end
end
