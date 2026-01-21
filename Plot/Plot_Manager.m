classdef Plot_Manager < handle
    %PLOT_MANAGER Centralized plotting style + export for a project
    
    properties
        % General
        x
        y_cell
        x_label_cell
        y_label_cell
        title_cell
        legend_cell

        print_legend
        filename
        save_pdf

        % Fixed parameters
        textwidth_cm = 13.75
        single_plot_position = [0.20, 0.17, 0.68, 0.72]

        num_Bins_hist = 30
        bar_width_hist = 0.8

        marker_size = 12

        scatter_size = 40
        scatter_linewidth = 1.3
        scatter_color = [0 0 0]

        fill_color = [0.8 0.8 1]
        fill_alpha = 0.35
    end
    
    methods
        function obj = Plot_Manager(args)
            %PLOT_MANAGER Construct an instance of this class
            % Assign values
            obj.x = args.x;
            obj.y_cell = args.y_cell;
            obj.x_label_cell = args.x_label_cell;
            obj.y_label_cell = args.y_label_cell;
            obj.title_cell = args.title_cell;
            obj.legend_cell = args.legend_cell;
            obj.print_legend = args.print_legend;
            obj.filename = args.filename;
            obj.save_pdf = args.save_pdf;
        end

        function axis_plot(obj, ax, i, opts)
            %render_axis creates plots per axis

            % Plot signals
            for j = 1:numel(obj.y_cell{i})
                plot(ax, obj.x, obj.y_cell{i}{j}, ...
                    'linewidth', opts.linewidth, ...
                    'Marker', opts.marker, ...
                    'MarkerSize', obj.marker_size);
            end
        end

        function axis_hist(obj, ax, x_hist)
            %axis_hist creates histogram in axis
            % Hist data
            [counts, edges] = histcounts(x_hist, obj.num_Bins_hist, 'Normalization', 'pdf');

            % Bin center
            binCenters = edges(1:end-1) + diff(edges)/2;

            % Plot hist
            bar(ax, binCenters, counts, ...
                'BarWidth', obj.bar_width_hist, ...
                'EdgeColor', 'none');
        end

        function axis_scatter(obj, ax, x_scatter, y_scatter)
            %axis_scatter creates scatter in axis
            % Plot scatter
            scatter(ax, x_scatter, y_scatter, obj.scatter_size, ...
                'Marker', '+', ...
                'LineWidth', obj.scatter_linewidth, ...
                'MarkerEdgeColor', obj.scatter_color)
        end

        function axis_fill(obj, ax, y_fill_upper, y_fill_lower)
            %axis_scatter creates fill area in axis
            % Fill
            hfill = fill(ax, [obj.x; flipud(obj.x)], [y_fill_upper; flipud(y_fill_lower)], ...
                obj.fill_color, ...
                'EdgeColor','none', ...
                'FaceAlpha', obj.fill_alpha);

            uistack(hfill, 'bottom')
        end
        
        function axis_options(obj, ax, i, opts)
            %axis_options set the axis options

            % General
            ax.Box = 'on';
            ax.XGrid = 'on';
            ax.YGrid = 'on';

            % Set y scale
            ax.YScale = opts.y_scale;

            % Set y limits
            yl = ax.YLim;
            dy = diff(yl);
            new_ylim = [yl(1) - opts.y_rel_offset*dy, yl(2) + opts.y_rel_offset*dy];
            ax.YLim = new_ylim;

            % Set x limits
            xl = [min(obj.x), max(obj.x)];
            dx = diff(xl);
            new_xlim = [xl(1) - opts.x_rel_offset*dx, xl(2) + opts.x_rel_offset*dx];
            ax.XLim = new_xlim;

            % Labels
            ax.XLabel.String = obj.x_label_cell{i};
            ax.XLabel.Interpreter = 'latex';
            ax.XLabel.FontSize = 11;
            
            ax.YLabel.String = obj.y_label_cell{i};
            ax.YLabel.Interpreter = 'latex';
            ax.YLabel.FontSize = 11;
            
            ax.Title.String = obj.title_cell{i};
            ax.Title.Interpreter = 'latex';
            ax.Title.FontSize = 12;

            % Legend
            if obj.print_legend
                legend(ax, obj.legend_cell{i}, 'Interpreter', 'latex', 'FontSize', 9, 'Location','best');
            end

            % Style Settings
            ax.Color = 'w';
            
            ax.FontName = 'Latin Modern Roman';
            ax.FontSize = 11;
            ax.TickLabelInterpreter = 'latex';

        end

        function export_plot(obj, fig, opts)
            %export_plot export plot to latex repo
            % Path to latex repo
            base_dir = fileparts(mfilename("fullpath"));
            Plot_base_path = fullfile(base_dir, '..', '..', 'Latex', 'Masterarbeit_Repo_Overleaf', 'Bilder', 'Kapitel');

            % Save Path
            Plot_path = fullfile(Plot_base_path, obj.filename);

            % Options
            fig.PaperUnits = 'centimeters';
            fig.PaperSize = [obj.textwidth_cm, opts.fig_height];
            fig.PaperPosition = [0 0 obj.textwidth_cm, opts.fig_height];

            if obj.save_pdf
                print(fig, Plot_path, '-dpdf');
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Plot Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function single_histo_plot(obj, opts, x_hist)
            %single_plot create plot with one axis in figure
            % Create figure
            fig = figure('Visible', 'on', ...
                           'Units', 'centimeters', ...
                           'Position', [2 2 obj.textwidth_cm opts.fig_height], ...
                           'Color', 'w');

            % Create axes inside figure
            ax = axes('Parent', fig);

            % Plot histogram
            hold(ax, 'on');
            obj.axis_hist(ax, x_hist);

            % Plot signals
            obj.axis_plot(ax, 1, opts);
            hold(ax, 'off');

            % Axis Options
            obj.axis_options(ax, 1, opts);

            % Axis size
            ax.Position = obj.single_plot_position;

            % Export figure
            obj.export_plot(fig, opts)
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
            hold(ax, 'on');
            obj.axis_plot(ax, 1, opts);
            hold(ax, 'off');

            % Axis Options
            obj.axis_options(ax, 1, opts);

            % Axis size
            ax.Position = obj.single_plot_position;

            % Export figure
            obj.export_plot(fig, opts)
        end

        function single_scatter_plot(obj, opts, x_scatter, y_scatter)
            %single_scatter_plot create plot with one axis in figure
            % Create figure
            fig = figure('Visible', 'on', ...
                           'Units', 'centimeters', ...
                           'Position', [2 2 obj.textwidth_cm opts.fig_height], ...
                           'Color', 'w');

            % Create axes inside figure
            ax = axes('Parent', fig);

            % Plot signals
            hold(ax, 'on');
            obj.axis_plot(ax, 1, opts);

            % Scatter
            obj.axis_scatter(ax, x_scatter, y_scatter);
            hold(ax, 'off');

            % Axis Options
            obj.axis_options(ax, 1, opts);

            % Axis size
            ax.Position = obj.single_plot_position;

            % Export figure
            obj.export_plot(fig, opts)
        end
        
        function single_scatter_fill_plot(obj, opts, x_scatter, y_scatter, y_fill_upper, y_fill_lower)
            %single_scatter_fill_plot create plot with one axis in figure
            % Create figure
            fig = figure('Visible', 'on', ...
                           'Units', 'centimeters', ...
                           'Position', [2 2 obj.textwidth_cm opts.fig_height], ...
                           'Color', 'w');

            % Create axes inside figure
            ax = axes('Parent', fig);

            % Plot signals
            hold(ax, 'on');
            obj.axis_plot(ax, 1, opts);

            % Scatter
            obj.axis_scatter(ax, x_scatter, y_scatter);

            % Fill area  
            obj.axis_fill(ax, y_fill_upper, y_fill_lower);
            hold(ax, 'off');

            % Axis Options
            obj.axis_options(ax, 1, opts);

            % Axis size
            ax.Position = obj.single_plot_position;

            % Export figure
            obj.export_plot(fig, opts)
        end

        function tiled_scatter_fill_plot(obj, opts, orientation, x_scatter_cell, y_scatter_cell, y_fill_upper_cell, y_fill_lower_cell)
            %tiled_scatter_fill_plot create plots with several axis in figure
            % Create figure
            fig = figure('Visible', 'on', ...
                           'Units', 'centimeters', ...
                           'Position', [2 2 obj.textwidth_cm opts.fig_height], ...
                           'Color', 'w');

            % Create tiled layout
            tiled = tiledlayout(fig, orientation(1), orientation(2), ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');

            % Create Plots in axes
            for i = 1:(orientation(1)*orientation(2))
                % Get ax object
                ax = nexttile(tiled);

                % Plot signals
                hold(ax, 'on');
                obj.axis_plot(ax, i, opts);
    
                % Scatter
                obj.axis_scatter(ax, x_scatter_cell{i}, y_scatter_cell{i});
    
                % Fill area  
                obj.axis_fill(ax, y_fill_upper_cell{i}, y_fill_lower_cell{i});
                hold(ax, 'off');
    
                % Axis Options
                obj.axis_options(ax, i, opts);
            end
        end
    end
end

