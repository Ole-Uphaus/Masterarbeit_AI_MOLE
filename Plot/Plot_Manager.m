classdef Plot_Manager < handle
    %PLOT_MANAGER Centralized plotting style + export for a project
    
    properties
        % General
        x
        y_cell
        x_label
        y_label_cell
        title_cell
        legend_cell

        filename
        save_pdf

        % Fixed parameters
        textwidth_cm = 13.75
    end
    
    methods
        function obj = Plot_Manager(args)
            %PLOT_MANAGER Construct an instance of this class
            % Assign values
            obj.x = args.x;
            obj.y_cell = args.y_cell;
            obj.x_label = args.x_label;
            obj.y_label_cell = args.y_label_cell;
            obj.title_cell = args.title_cell;
            obj.legend_cell = args.legend_cell;
            obj.filename = args.filename;
            obj.save_pdf = args.save_pdf;
        end

        function axis_plot(obj, ax, i, opts)
            %render_axis creates plots per axis

            % Plot signals
            hold(ax, 'on');
            for j = 1:numel(obj.y_cell{i})
                plot(ax, obj.x, obj.y_cell{i}{k}, 'linewidth', opts.linewidth);
            end
            hold(ax, 'off');
        end
        
        function axis_options(obj, ax, opts)
            %axis_options set the axis options

            % General
            ax.Box = 'on';
            ax.XGrid = 'on';
            ax.YGrid = 'on';

            % Set limits
            yl = ax.YLim;
            dy = diff(yl);
            new_ylim = [yl(1) - 0.05*dy, yl(2) + 0.05*dy];
            ax.YLim = new_ylim;

            % Style Settings
            ax.Color = 'w';
            
            ax.FontName = 'Latin Modern Roman';
            ax.FontSize = 11;
            ax.TickLabelInterpreter = 'latex';

            ax.Position = [0.20, 0.17, 0.68, 0.72];

        end

        function single_plot(obj, opts)
            %single_plot create plot with one axis in figure
            % Create figure
            fig = figure('Visible', 'on', ...
                           'Units', 'centimeters', ...
                           'Position', [2 2 obj.textwidth_cm opts.fig_height], ...
                           'Color', 'w');

            % Create axes inside figure
            ax = axes('Parent', fig);

            % Plot signals
            obj.axis_plot(ax, 1, opts);


        end      
    end
end

