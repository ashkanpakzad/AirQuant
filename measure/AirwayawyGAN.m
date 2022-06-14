classdef AirwayawyGAN < SuperMeasure
    methods
        function obj = AirwayawyGAN(varargin)
            % See superclass constructor method SuperMeasure.
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
            obj@SuperMeasure(varargin{:});
        end

        function Measure(obj, modulepath, model_path)
            % Measure airways using a deep learning method
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   modulepath(char): absolute path to awyGAN parent directory.
            %   model_path(char): absolute path to model load checkpoint.
            %
            %
            

            % save patches to temp dir
            obj.pixsize = 0.5;
            temp_dir = tempname;
            obj.tube.ExportOrthoPatches(temp_dir, 'temp')

            % compute image stats
            M = cell2mat(obj.tube.source);
            meanval = mean(M(:));
            varval = std(M(:));
            
            % add module path to python paths and import
            GANCNNpath = fullfile(modulepath,'GANCNN');
            pyrun(['import sys; sys.path.append("', modulepath,'");']);
            pyrun(['import sys; sys.path.append("', GANCNNpath,'");']);
            mod = py.importlib.import_module('AQ_CNR');
            py.importlib.reload(mod);

            pycode = ...
                ['import sys; sys.path.append("', modulepath,'");', ...
                ' from AQ_CNR import runawyGANmodel;',...
                ' modelpath = "', model_path, '";', ...
                ' data_path = "', temp_dir, '";', ...
                ' batch_size = ', num2str(256), ';', ...
                ' mean = ', num2str(meanval), ';', ...
                ' var = ', num2str(varval), ';', ...
                ' out = runawyGANmodel(modelpath, data_path, batch_size, mean, var)'];
            
            pymeasures = pyrun(pycode, 'out');
            measures = double(pymeasures);

            % set up memory
            slices_sz = length(obj.tube.source);
            ellipse_inner = cell(slices_sz, 1);
            ellipse_outer = cell(slices_sz, 1);
            
            for k = 1:slices_sz
                imcc = size(obj.tube.source{k}) / 2 - obj.pixsize;
                center = imcc + measures(k,3:4) / obj.pixsize;

                % process inner ellipses
                innerstruct = struct;
                innerstruct.center = center;
                innerstruct.Rx = measures(k,1)/obj.pixsize(1);
                innerstruct.Ry = measures(k,2)/obj.pixsize(1);
                innerstruct.rotation = measures(k,7);

                % process outer ellipses
                outerstruct = struct;
                outerstruct.center = center;
                outerstruct.Rx = measures(k,5)/obj.pixsize(1);
                outerstruct.Ry = measures(k,6)/obj.pixsize(1);
                outerstruct.rotation = measures(k,7);

                % compute ellipses
                ellipse_inner{k,1} = AQEllipse(obj.pixsize, innerstruct);
                ellipse_outer{k,1} = AQEllipse(obj.pixsize, outerstruct);
            end
            obj.measures = cell(2, length(ellipse_inner));
            obj.measures(1,:) = ellipse_inner;
            obj.measures(2,:) = ellipse_outer;
            
        end

        
    end
end
