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
            % view in itksnap
            %
            % May need to set up enviroments on matlab search path for
            % system terminal.
            % setenv('PATH', [getenv('PATH') ':/Applications/ITK-SNAP.app/Contents/bin']);
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
            % short desc
            %
            % long desc
            %
            % .. todo::
            %   * add documentation to this function
            %   * add version that makes tubestack back to parent
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %

            arguments
                obj
                options.type
            end
            volout = get(obj,options.type);
        end
    end
    methods(Static)
        function list = list_property(someobject, property)
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