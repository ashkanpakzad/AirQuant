% Lead Author: Ashkan Pakzad 2021. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef AirQuant < handle & matlab.mixin.SetGet & matlab.mixin.Copyable
    properties
        source
        seg
        voxdim
    end
    methods
        function obj = AirQuant()
        end
        function toNii(obj, filename, options)
            % Save 3D array object property as nifti image.
            %
            % Args:
            %   filename(`char`): name to save as.
            %   type(char): *OPTIONAL* `default = 'source'` object
            %       property that is a 3D array.
            %   gz(bool): *OPTIONAL* `default = true` whether to save gzip
            %       compress the output.
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> AQnet.toNii('example');
            %

            arguments
                obj
                filename char
                options.type = 'source'
                options.gz logical = true
            end

            filename = parse_filename_extension(filename, '.nii');

            volout = ParseVolOut(obj, type=options.type);
            % nifti convention typically limited to int16
            volout = int16(volout);
            if isprop(obj,'sourceinfo')
                info = obj.sourceinfo;
            elseif isprop(obj,'network')
                info = obj.network.sourceinfo;
            end
            info.ImageSize = size(volout);
            info.PixelDimensions = obj.voxdim;
            info.Datatype = 'int16';
            info.Description = 'Output from AirQuant';

            niftiwrite(volout,filename,info)

            if options.gz == true
                gzip(filename)
                delete(filename)
            end
        end

        function status = toITKsnap(obj, segname)
            % View images in `ITK-snap <http://www.itksnap.org>`_ medical image viewer.
            %
            % Open source image and segmentation overlay in 
            % `ITK-snap (www.itksnap.org) <http://www.itksnap.org>`_.
            % Images are saved as temporary files and then called using
            % system commands to open in ITK-snap.
            %
            % .. note:
            %   `ITK-snap (www.itksnap.org) <http://www.itksnap.org>`_
            %   needs to be installed on the system and on both the system
            %   and matlab search path. e.g. ``` setenv('PATH',
            %   [getenv('PATH')
            %   ':/Applications/ITK-SNAP.app/Contents/bin']);```
            %
            % Args:
            %   type(char): *OPTIONAL* `default = 'seg'`. object
            %       property name of a segmentation.
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> AQnet.toITKsnap();
            %

            if nargin < 2
                segname = 'seg';
            end

            % set up temp file
            sourcefile = parse_filename_extension(tempname, '.nii');
            obj.toNii(sourcefile, type='source', gz=false)
            command = ['itksnap ', sourcefile];

            if ~isempty(segname)
                segfile = parse_filename_extension(tempname, '.nii');
                obj.toNii(segfile, type=segname, gz=false)
                command = ['itksnap -g ', sourcefile, ' -s ', segfile];
            end

            status = system(command);

            if status ~= 0
                error(['Failed to open in ITK-snap, ' ...
                    'please see documentation to check this facility is ' ...
                    'set up correctly.'])
            end

        end

        function volout = ParseVolOut(obj,options)
            % Parse volume of a given property name
            %
            % Parses the volume of a given property name as output.
            %
            %
            % Args:
            %   type(char): *OPTIONAL* `default = 'source'`. Object
            %       property name.
            %
            % Return:
            %   1 variable
            %   * volout(`3D array`): output volume.
            %
            % .. todo:
            %   Needs example.
            %
            % Example:
            %   >>> run CA_base.m;
            %

            arguments
                obj
                options.type = 'source'
            end
            volout = get(obj,options.type);
        end
    end
    methods(Static)
        function list = list_property(someobject, property)
            % Get the same property for all objects in a list.
            %
            % Make an array of an object's property given an array of that
            % object.
            %
            %
            % Args:
            %   someobject(object array): Object cell array.
            %
            % Return:
            %   1 variable
            %   * list(`array`): output object property array, where '' is
            %       used if that property for a particular object does not exist.
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> sourcevol = AQnet.ParseVolOut(type='source');
            list = cellfun(@nanifempty, someobject, 'UniformOutput', false);

            function out = nanifempty(c)
                if isfield(c,property)
                    out = c.(property);
                else
                    out = '' ;
                end
            end
        end
    end
end