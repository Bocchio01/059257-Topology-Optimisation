classdef Node
    properties
        x
        y
        forceX = 0
        forceY = 0
        isConstrainedX = false
        isConstrainedY = false
    end
    
    methods
        function obj = Node(x, y, isConstrainedX, isConstrainedY, forceX, forceY)
            % Constructor for Node class
            obj.x = x;
            obj.y = y;
            
            if nargin > 2 % if constraints are provided
                obj.isConstrainedX = isConstrainedX;
                obj.isConstrainedY = isConstrainedY;
            end
            
            if nargin > 4 % if forces are provided
                obj.forceX = forceX;
                obj.forceY = forceY;
            end
        end
        
        function coords = getCoordinates(obj)
            % Returns the coordinates of the node
            coords = [obj.x, obj.y];
        end
        
        function force = getForce(obj)
            % Returns the force vector at the node
            force = [obj.forceX, obj.forceY];
        end
        
        function isEqual = eq(obj1, obj2)
            % Custom equality check for Node objects

            isEqual = ([obj1.x] == obj2.x) & ([obj1.y] == obj2.y);
        end
    end
end
