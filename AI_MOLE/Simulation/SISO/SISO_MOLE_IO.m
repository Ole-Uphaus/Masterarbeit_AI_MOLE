classdef SISO_MOLE_IO < handle
    %SISO_MOLE_IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = SISO_MOLE_IO()
            %SISO_MOLE_IO Construct an instance of this class
            %   Detailed explanation goes here
            
            % Geerate dynamic paths
            base_dir = fileparts(mfilename("fullpath"));
            GP_path = fullfile(base_dir, '..', '..', '..', 'GP', 'GP_SISO');
            ILC_path = fullfile(base_dir, '..', '..', '..', 'ILC', 'ILC_SISO');

            % Add temporary paths
            addpath(GP_path);
            addpath(ILC_path);
        end

    end
end

