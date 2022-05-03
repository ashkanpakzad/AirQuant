classdef HumanAirways < TubeNetwork
    % short desc
    %
    % long desc
    %
    % .. todo: add documentation to this function
    %
    % Args:
    %   x(type):
    %
    % Return:
    %   y(type):
    %

    properties
        Property1
    end

    methods
        function obj = HumanAirways(inputArg1,inputArg2)
            % short desc
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = ParseSeg(obj,inputArg)
            % short desc
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            outputArg = obj.Property1 + inputArg;
        end
        
        function outputArg = MakeDigraph(obj,inputArg)
            % short desc
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            outputArg = obj.Property1 + inputArg;
        end
        
        function IdentifyCarinaAndTrachea(obj)
            % short desc
            %
            % long desc
            %
            % .. todo:
            %   * May need to add further methods to reidentify generations
            %   for all airways.
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %

            g = TubesAsEdges(obj);
            [~, carina_node] = max(centrality(g,'outcloseness'));
            carinaend_gidx = inedges(carina_node);
            obj.tubes(g.Edges.ID(carinaend_gidx)).carinaend = true;
            tracheanodes = predecessors(g, carina_node);
            for nid = tracheanodes
                eid = inedges(g, nid);
                tubeid = g.Edges(eid).ID;
                obj.tubes(tubeid).istrachea = true;
                obj.tubes(tubeid).generation = 0;
            end
        end

        function outputArg = IdentifyTrachea(obj,inputArg)
            % short desc
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            outputArg = obj.Property1 + inputArg;
        end
        
        function outputArg = ClassifyLobes(obj,inputArg)
            % short desc
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            outputArg = obj.Property1 + inputArg;
        end
    end
end

