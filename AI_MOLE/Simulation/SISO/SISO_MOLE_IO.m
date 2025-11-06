classdef SISO_MOLE_IO < handle
    %SISO_MOLE_IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        GP_SISO
        ILC_SISO
        y_cell
        u_cell
    end
    
    methods
        function obj = SISO_MOLE_IO(r_vec, m_delay, u_init)
            %SISO_MOLE_IO Construct an instance of this class
            %   Detailed explanation goes here
            
            % Generate dynamic paths
            base_dir = fileparts(mfilename("fullpath"));
            GP_path = fullfile(base_dir, '..', '..', '..', 'GP', 'GP_SISO');
            ILC_path = fullfile(base_dir, '..', '..', '..', 'ILC', 'ILC_SISO');

            % Add temporary paths
            addpath(GP_path);
            addpath(ILC_path);

            % Initialize Classes
            obj.GP_SISO = GP_SISO_IO();
            obj.ILC_SISO = ILC_SISO(r_vec, m_delay, u_init);
        end

        function update_input(obj, y_vec)
        end
    end
end

