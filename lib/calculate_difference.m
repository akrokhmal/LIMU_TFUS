function metrics = calculate_difference(field1_name, field1, field2_name, field2, dx, z_start, varargin)
%COMPUTEDIFFERENCEMETRICS Compute difference between two fields
%
% DESCRIPTION:
%     computeDifferenceMetrics computes various difference metrics to
%     quantify the difference between two scalar fields. For calculating
%     relative differences, field 1 is always taken as the reference field.
%
%     The input fields can be either 2D or 3D. If the fields are 2D, the
%     difference in focal volume is calculated as a % difference in area.
%     If 3D, it is calculated as a % difference in volume. 
%
%     To exclude the first part of the field from the comparison, set the
%     value for exit_plane_ind. 
%
%     The error norms are calculated over the entire field (ignoring the
%     mask), with percentage values taken relative to the maximum value in
%     field 1. The focal differences are calculated over the field after
%     accounting for the mask. 1D profiles are taken through the maximum
%     value in the masked region of field 1.
%
%     The fields and differences can optionally be displayed by setting
%     plot_fields to true. The names are used for plotting only.
%
% USAGE:
%     metrics = computeDifferenceMetrics(field1, field2, plot_fields)
%
% INPUTS:
%     field1_name       - Name of first field.
%     field1            - Values of first field.
%     field2_name       - Name of second field.
%     field2            - Values Second field.
%     dx                - Grid spacing. [m]
%     z_start - starting z coordinate [m]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'ExitPlaneIndex'  - Index of the exit plane in the x-dim. Used to
%                         exclude part of the field from the metrics
%                         (default = 1).
%     'LateralProfileIndex'
%                       - Sets the position of the maximum in the x
%                         direction. Useful for comparing fields in which
%                         there is not a natural focus.
%     'FocalMask'       - Mask used for calculation of the peak pressure in
%                         the region of interest (default = complete
%                         field).
%     'Plot'            - Boolean controlling whether fields are plotted
%                         (default = true).
%     'VolumeMetrics'   - Boolean controlling whether volume metrics are
%                         computed (default = false). Requires the Image
%                         Processing Toolbox.
%
% OUTPUTS:
%     metrics           - Structure containing the following fields:
%
%                         linf_perc:             relative linf error [%]
%                         l2_perc:               relative l2 error [%]
%                         max_amp_perc:          difference in max amplitude [%]
%                         max_pos_mm:            difference in max position [m]
%                         focal_size_3dB_x_mm:   difference in focal size in x-direction [m]
%                         focal_size_3dB_y_mm:   difference in focal size in y-direction [m]
%                         focal_size_3dB_z_mm:   difference in focal size in z-direction [m]
%                         focal_size_6dB_x_mm:   difference in focal size in x-direction [m]
%                         focal_size_6dB_y_mm:   difference in focal size in y-direction [m]
%                         focal_size_6dB_z_mm:   difference in focal size in z-direction [m]
%
%                         If 'VolumeMetrics', true:
%
%                         focal_volume_6dB_mm2:    difference in focal volume (2D)
%                         focal_volume_6dB_mm3:    difference in focal volume (3D)
%                         focal_volume_6dB_perc:   difference in focal volume
%
% ABOUT:
%     author            - Bradley Treeby
%     date              - 27 February 2021
%     last update       - Alisa Krokhmal, 03 Sep 2024

% set default values
num_req_input_variables = 5;
exit_plane_ind          = 1;
mask                    = ones(size(field1));
plot_fields             = true;
lateral_profile_x_ind   = false;
compute_volume_metrics  = false;

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'ExitPlaneIndex'
                exit_plane_ind = varargin{input_index + 1}; 
            case 'FocalMask'
                mask = varargin{input_index + 1}; 
            case 'LateralProfileIndex'
                lateral_profile_x_ind = varargin{input_index + 1}; 
            case 'Plot'
                plot_fields = varargin{input_index + 1};
            case 'VolumeMetrics'
                compute_volume_metrics = varargin{input_index + 1};                
            otherwise
                error(['Unknown optional input: ' varargin{input_index} '.']);
        end
    end
end

% get grid size
[Nx, Ny, Nz] = size(field1);

% account for exit plane
field1(1:exit_plane_ind - 1, :, :) = 0;
field2(1:exit_plane_ind - 1, :, :) = 0;

% maximum linf over the whole field
diff = 100 * abs(field1 - field2) ./ max(abs(field1(:)));
metrics.linf_perc = max(diff(:));

% relative l2 over the whole filed 
metrics.l2_perc = computeL2(field1, field2);

metrics.p1_max = maxND(field1);
metrics.p2_max = maxND(field2);

% get maximum value and the position of the maximum within the masked area
field1_masked = field1;
field1_masked(mask == 0) = 0;
field2_masked = field2;
field2_masked(mask == 0) = 0;
if lateral_profile_x_ind
    [mx1, ind_yz1] = maxND(abs(field1_masked(lateral_profile_x_ind, :, :)));
    [mx2, ind_yz2] = maxND(abs(field2_masked(lateral_profile_x_ind, :, :)));
    ps1 = [lateral_profile_x_ind, ind_yz1(2:end)];
    ps2 = [lateral_profile_x_ind, ind_yz2(2:end)];
else

    [mx1, ps1] = maxND(abs(field1_masked));
    [mx2, ps2] = maxND(abs(field2_masked));
end
if Nz == 1
    ps1 = [ps1, 1];
    ps2 = [ps2, 1];
end

% difference in peak pressure within the masked area
metrics.max_amp_perc = 100 * abs(mx1 - mx2) / abs(mx1);

% difference in the position of the peak
metrics.max_pos_mm = 1e3 * dx * norm(ps1 - ps2, 2);

metrics.max_pos1 = 1e3 * dx * ps1;
metrics.max_pos1(1) = metrics.max_pos1(1) - 1e3*dx*(Nx-1)/2;
metrics.max_pos1(2) = metrics.max_pos1(2) - 1e3*dx*(Ny-1)/2;
metrics.max_pos1(3) = metrics.max_pos1(3) + 1e3*z_start;
metrics.max_pos2 = 1e3 * dx * ps2;
metrics.max_pos2(1) = metrics.max_pos2(1) - 1e3*dx*(Nx-1)/2;
metrics.max_pos2(2) = metrics.max_pos2(2) - 1e3*dx*(Ny-1)/2;
metrics.max_pos2(3) = metrics.max_pos2(3) + 1e3*z_start;


% axial and lateral profiles through the peak (always using the
% position of the peak for the reference field)
field1_profile_x = field1(:, ps1(2), ps1(3));
field2_profile_x = field2(:, ps1(2), ps1(3));
field1_profile_y = field1(ps1(1), :, ps1(3));
field2_profile_y = field2(ps1(1), :, ps1(3));
field1_profile_z = field1(ps1(1), ps1(2), :);
field2_profile_z = field2(ps1(1), ps1(2), :);

% -3 dB focal size
plot_fwhm = false;
field1_focal_size_3dB_x = fwhm2(field1_profile_x, dx, ps1(1), plot_fwhm);
field2_focal_size_3dB_x = fwhm2(field2_profile_x, dx, ps2(1), plot_fwhm);
field1_focal_size_3dB_y = fwhm2(field1_profile_y, dx, ps1(2), plot_fwhm);
field2_focal_size_3dB_y = fwhm2(field2_profile_y, dx, ps2(2), plot_fwhm);
field1_focal_size_3dB_z = fwhm2(field1_profile_z, dx, ps1(3), plot_fwhm);
field2_focal_size_3dB_z = fwhm2(field2_profile_z, dx, ps2(3), plot_fwhm);

% -6 dB focal size
field1_focal_size_6dB_x = fwhm2(field1_profile_x.^2, dx, ps1(1), plot_fwhm);
field2_focal_size_6dB_x = fwhm2(field2_profile_x.^2, dx, ps2(1), plot_fwhm);
field1_focal_size_6dB_y = fwhm2(field1_profile_y.^2, dx, ps1(2), plot_fwhm);
field2_focal_size_6dB_y = fwhm2(field2_profile_y.^2, dx, ps2(2), plot_fwhm);
field1_focal_size_6dB_z = fwhm2(field1_profile_z.^2, dx, ps1(3), plot_fwhm);
field2_focal_size_6dB_z =fwhm2(field2_profile_z.^2, dx, ps2(3), plot_fwhm);

% difference in focal size
metrics.focal_size_3dB_x_mm = 1e3 * (field1_focal_size_3dB_x - field2_focal_size_3dB_x);
metrics.focal_size_3dB_y_mm = 1e3 * (field1_focal_size_3dB_y - field2_focal_size_3dB_y);
metrics.focal_size_3dB_z_mm = 1e3 * (field1_focal_size_3dB_z - field2_focal_size_3dB_z);

metrics.focal_size_6dB_x_mm = (field1_focal_size_6dB_x - field2_focal_size_6dB_x);
metrics.focal_size_6dB_y_mm = (field1_focal_size_6dB_y - field2_focal_size_6dB_y);
metrics.focal_size_6dB_z_mm = (field1_focal_size_6dB_z - field2_focal_size_6dB_z);

metrics.focal_zone1 = [field1_focal_size_6dB_x, field1_focal_size_6dB_y, field1_focal_size_6dB_z]*1e3;
metrics.focal_zone2 = [field2_focal_size_6dB_x, field2_focal_size_6dB_y, field2_focal_size_6dB_z]*1e3;

% calculations for z direction if 3D
if Nz ~= 1
    field1_profile_z = squeeze(field1(ps1(1), ps1(2), :));
    field2_profile_z = squeeze(field2(ps1(1), ps1(2), :));
    
    field1_focal_size_3dB_z = fwhm2(field1_profile_z, dx, ps1(3), plot_fwhm);
    field2_focal_size_3dB_z = fwhm2(field2_profile_z, dx, ps2(3), plot_fwhm);
    
    field1_focal_size_6dB_z = fwhm2(field1_profile_z.^2, dx, ps1(3), plot_fwhm);
    field2_focal_size_6dB_z = fwhm2(field2_profile_z.^2, dx, ps2(3), plot_fwhm);
   
    metrics.focal_size_3dB_z_mm = 1e3 * (field1_focal_size_3dB_z - field2_focal_size_3dB_z);
    
    metrics.focal_size_6dB_z_mm = 1e3 * (field1_focal_size_6dB_z - field2_focal_size_6dB_z);
else
    metrics.focal_size_3dB_z_mm = 0;
    metrics.focal_size_6dB_z_mm = 0;
end

if compute_volume_metrics

    % -6dB focal volume for field 1
    CC = bwconncomp(field1_masked > 0.5 * max(field1_masked(:)));
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [num_vox_focal_volume1, idx] = max(numPixels);
    focal_volume1 = zeros(size(field1));
    focal_volume1(CC.PixelIdxList{idx}) = 1;

    % -6dB focal volume for field 2
    CC = bwconncomp(field2_masked > 0.5 * max(field2_masked(:)));
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [num_vox_focal_volume2, idx] = max(numPixels);
    focal_volume2 = zeros(size(field1));
    focal_volume2(CC.PixelIdxList{idx}) = 1;

    % difference in -6dB focal volume
    focal_volume_6dB   = (num_vox_focal_volume1 - num_vox_focal_volume2);
    metrics.focal_volume_6dB_perc = 100 * abs(num_vox_focal_volume2 - num_vox_focal_volume1) / num_vox_focal_volume1;
    if Nz == 1
        metrics.focal_volume_6dB_mm2 = 1e6 * focal_volume_6dB * dx^2;
    else
        metrics.focal_volume_6dB_mm3 = 1e9 * focal_volume_6dB * dx^3;
    end
    
end

metrics.vol1 = dx^3*1e9*num_vox_focal_volume1;
metrics.vol2 = dx^3*1e9*num_vox_focal_volume2;


if plot_fields

    % define axis vectors for plotting %%%% CHANGED BY ALISA
    z_vec = 1e3*z_start + 1e3 * (0:(Nz-1)) * dx;
    y_vec = 1e3 * (-(Ny-1)/2:(Ny-1)/2) * dx;
    x_vec = 1e3 * (-(Nx-1)/2:(Nx-1)/2) * dx;
    
    % color scale
    c_sc = [min(field1(:)), max(field1(:))] * 1e-3;
    
    % plot fields (planes through peak of reference field in 3D)
    for ind = 1:3
            
        switch ind
            case 1
                field1_plot = field1(:, :, ps1(3));
                field2_plot = field2(:, :, ps1(3));
                plot_dim1 = 'X';
                plot_dim2 = 'Y';
                ax1 = x_vec;
                ax2 = y_vec;
            case 2
                field1_plot = squeeze(field1(:, ps1(2), :));
                field2_plot = squeeze(field2(:, ps1(2), :));
                plot_dim1 = 'X';
                plot_dim2 = 'Z';
                ax1 = x_vec;
                ax2 = z_vec;                    
            case 3
                field1_plot = squeeze(field1(ps1(1), :, :));
                field2_plot = squeeze(field2(ps1(1), :, :));
                plot_dim1 = 'Y';
                plot_dim2 = 'Z';
                ax1 = y_vec;
                ax2 = z_vec;
        end

        figure;
        subplot(1, 3, 1);
        imagesc(ax2, ax1, 1e-3 * field1_plot, c_sc);
        xlabel(['Lateral Position (' plot_dim2 ') [mm]']);
        ylabel(['Axial Position (' plot_dim1 ') [mm]']);
        axis image;
        cb = colorbar;
        title(cb, '[kPa]');
        title([field1_name ' (REFERENCE)']);

        subplot(1, 3, 2);
        imagesc(ax2, ax1, 1e-3 * field2_plot, c_sc);
        xlabel(['Lateral Position (' plot_dim2 ') [mm]']);
        ylabel(['Axial Position (' plot_dim1 ') [mm]']);
        axis image;
        cb = colorbar;
        title(cb, '[kPa]');
        title(field2_name);

        subplot(1, 3, 3);
        imagesc(ax2, ax1, 100 * abs(field1_plot - field2_plot) / max(field1_plot(:)));
        xlabel(['Lateral Position (' plot_dim2 ') [mm]']);
        ylabel(['Axial Position (' plot_dim1 ') [mm]']);
        axis image;
        colorbar;
        cb = colorbar;
        title(cb, '[%]');
        title('DIFFERENCE');

        sgtitle([field2_name ' - Pressure Amplitude (' plot_dim1 plot_dim2 ')']);
        scaleFig(2, 1);
        drawnow;

        if compute_volume_metrics
            
            switch ind
                case 1
                    volume1_plot = focal_volume1(:, :, ps1(3));
                    volume2_plot = focal_volume2(:, :, ps1(3));
                case 2
                    volume1_plot = squeeze(focal_volume1(:, ps1(2), :));
                    volume2_plot = squeeze(focal_volume2(:, ps1(2), :));
                case 3
                    volume1_plot = squeeze(focal_volume1(ps1(1), :, :));
                    volume2_plot = squeeze(focal_volume2(ps1(1), :, :));
            end

            figure;
            subplot(1, 3, 1);
            imagesc(ax2, ax1, volume1_plot);
            xlabel(['Lateral Position (' plot_dim2 ') [mm]']);
            ylabel(['Axial Position (' plot_dim1 ') [mm]']);
            axis image;
            cb = colorbar;
            title(cb, '[kPa]');
            title([field1_name ' (REFERENCE)']);

            subplot(1, 3, 2);
            imagesc(ax2, ax1, volume2_plot);
            xlabel(['Lateral Position (' plot_dim2 ') [mm]']);
            ylabel(['Axial Position (' plot_dim1 ') [mm]']);
            axis image;
            cb = colorbar;
            title(cb, '[kPa]');
            title(field2_name);

            subplot(1, 3, 3);
            imagesc(ax2, ax1, abs(volume1_plot - volume2_plot));
            xlabel(['Lateral Position (' plot_dim2 ') [mm]']);
            ylabel(['Axial Position (' plot_dim1 ') [mm]']);
            axis image;
            colorbar;
            cb = colorbar;
            title(cb, '[%]');
            title('DIFFERENCE');

            sgtitle([field2_name ' - Focal Volume (' plot_dim1 plot_dim2 ')']);
            scaleFig(2, 1);
            drawnow;
            
        end

        if Nz == 1
            break;
        end
        
    end
        
    % plot profiles through peak of reference field
    if Nz == 1
        num_cols = 2;
    else
        num_cols = 3;
    end
    
    figure;
    subplot(2, num_cols, 1);
    plot(x_vec, 1e-3 * field1_profile_x, 'LineWidth', 1.5);
    hold on;
    plot(x_vec, 1e-3 * field2_profile_x, '-', 'LineWidth', 1.5);
    grid on;
    legend(field1_name, field2_name, 'Location', 'Best');
    xlabel('Axial Position (X) [mm]');
    ylabel('Pressure [kPa]');
    title(['X Profile (Y = ' num2str(y_vec(ps1(2))) ' mm, Z = ' num2str(z_vec(ps1(3))) ' mm)']);
    
    subplot(2, num_cols, 2);
    plot(y_vec, 1e-3 * field1_profile_y, 'LineWidth', 1.5);
    hold on;
    plot(y_vec, 1e-3 * field2_profile_y, '-', 'LineWidth', 1.5);
    grid on;
    xlabel('Lateral Position (Y) [mm]');
    ylabel('Pressure [kPa]');
    title(['Y Profile (X = ' num2str(x_vec(ps1(1))) ' mm, Z = ' num2str(z_vec(ps1(3))) ' mm)']);
    
    if Nz > 1
        subplot(2, num_cols, 3);
        plot(z_vec, 1e-3 * field1_profile_z, 'LineWidth', 1.5);
        hold on
        plot(z_vec, 1e-3 * field2_profile_z, '-', 'LineWidth', 1.5);
        grid on;
        xlabel('Lateral Position (Z) [mm]');
        ylabel('Pressure [kPa]');
        title(['Z Profile (X = ' num2str(x_vec(ps1(1))) ' mm, Y = ' num2str(y_vec(ps1(2))) ' mm)']);
    end
    
    subplot(2, num_cols, num_cols + 1);
    plot(x_vec, 100 * abs(field1_profile_x - field2_profile_x) ./ max(field1_profile_x(:)), 'LineWidth', 1.5);
    grid on;
    xlabel('Axial Position (X) [mm]');
    ylabel('Difference [% of profile peak]');
    title('Axial Profile (X) Difference');
    
    subplot(2, num_cols, num_cols + 2);
    plot(y_vec, 100 * abs(field1_profile_y - field2_profile_y) ./ max(field1_profile_y(:)), 'LineWidth', 1.5);
    grid on;
    xlabel('Lateral Position (Y) [mm]');
    ylabel('Difference [% of profile peak]');
    title('Lateral Profile (Y) Difference');
    
    if Nz > 1
        subplot(2, num_cols, num_cols + 3);
        plot(z_vec, 100 * abs(field1_profile_z - field2_profile_z) ./ max(field1_profile_z(:)), 'LineWidth', 1.5);
        grid on;
        xlabel('Lateral Position (Y) [mm]');
        ylabel('Difference [% of profile peak]');
        title('Lateral Profile (Z) Difference');
    end
    
    vol_diff = '';
    if Nz > 1
        fwhm_diff = ['\Delta FWHM = (' num2str(metrics.focal_size_6dB_x_mm, 2) ', ' ...
            num2str(metrics.focal_size_6dB_y_mm, 2) ', ' ...
            num2str(metrics.focal_size_6dB_z_mm, 2) ') mm'];
        
        if compute_volume_metrics
            vol_diff = [' || \Delta V = ' num2str(metrics.focal_volume_6dB_mm3, 2) ' mm^3 (' num2str(metrics.focal_volume_6dB_perc, 2) '%)'];
        end
    else
        fwhm_diff = ['\Delta FWHM = (' num2str(metrics.focal_size_6dB_x_mm, 2) ', ' ...
            num2str(metrics.focal_size_6dB_y_mm, 2) ') mm'];
        
        if compute_volume_metrics
            vol_diff = [' || \Delta A = ' num2str(metrics.focal_volume_6dB_mm2, 2) ' mm^2 (' num2str(metrics.focal_volume_6dB_perc, 2)  '%)'];
        end
    end
    
    sgtitle({field2_name, [...
    '\Delta p_{focus} = '   num2str(metrics.max_amp_perc, 2) '% || ' ...
    '\Delta focus pos = '   num2str(metrics.max_pos_mm, 2) ' mm || ' ...
    'L_\infty = '           num2str(metrics.linf_perc, 2) '% || ', ...
    'L2 = '                num2str(metrics.l2_perc, 2) '% || ', ...
    fwhm_diff, vol_diff]});
    
    scaleFig(2, 1);
    drawnow;
    
end

% -------------

function l2 = computeL2(field1, field2)
% Function to compute normalised L2 error ignoring nans in either field.

not_nan_index = ~(isnan(field1) | isnan(field2));
l2 = 100 * sqrt( sum((field1(not_nan_index) - field2(not_nan_index)).^2) ./ sum(field1(not_nan_index).^2) );

% -------------

function fwhm_val = fwhm2(f, x, i_max, plot_fwhm)
%FWHM Compute the full width at half maximum.
%
% DESCRIPTION:
%     fwhm calculates the Full Width at Half Maximum (FWHM) of a positive
%     1D input function f(x) with spacing given by x.
%
% USAGE:
%     fwhm_val = fwhm2(f, x, i_max, plot_fwhm)
%
% INPUTS:
%     f           - f(x) data
%     x           - x data or dx
%     i_max       - grid index of focus
%
% OPTIONAL INPUTS:
%     plot_fwhm   - Boolean controlling whether a plot of the function
%                   and FWHM is produced.
%
% OUTPUTS:
%     fwhm_val    - FWHM of f(x)
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 5th May 2009
%     last update - 2nd July 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2021 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

% check for optional inputs
if nargin == 3
    plot_fwhm = false;
end

% check if dx is given in place of an x array
if length(x) == 1
    x = 0:x:(length(f) - 1) * x;
end

% check the input is 1D
if numDim(f) ~= 1
    error('Input function must be 1 dimensional.');
end

% find the maximum value of f(x)
f_max = f(i_max);

% setup the indexing variables for finding the leading edge
index = i_max;

% loop until the index at half maximum is found
while f(index) > (0.5 * f_max)
    index = index - 1;
    if index < 1
        warning('Left half maximum not found.');
        fwhm_val = nan;
        return
    end
end

% linearly interpolate between the previous values to find the value of x
% at the leading edge at half maximum
m = (f(index+1) - f(index)) / (x(index+1) - x(index));
c = f(index) - x(index) * m;
i_leading = (0.5 * f_max  - c) / m;

% setup the indexing variables for finding the trailing edge
index = i_max;

% loop until the index at half maximum is found
while f(index) > (0.5 * f_max)
    index = index + 1;
    if index > length(f)
        warning('Right half maximum not found.');
        fwhm_val = nan;
        return
    end
end

% linearly interpolate between the previous values to find the value of x
% at the trailing edge at half maximum
m = (f(index - 1) - f(index)) / (x(index - 1) - x(index));
c = f(index) - x(index) * m;
i_trailing = (0.5 * f_max  - c) / m;
    
% compute the FWHM
fwhm_val = abs(i_trailing - i_leading);

% plot the function and the FWHM if required
if plot_fwhm   
    figure;
    plot(x, f, 'b-');
    xlabel('x');
    ylabel('f(x)');    
    hold on;
    plot([i_leading, i_trailing], [0.5 * f_max, 0.5 * f_max], 'r*-');
    title(['Full Width At Half Maximum (FWHM): ' num2str(fwhm_val) ]);
end
