function vHolo = bowl_source_hologram_calculation(c0, rho0, a, r, f, dx, v0, source_pos, varargin)

%DESCRIPTION:
% bowl_source_hologram_calculation computes the hologram of an ideal 
% bowl- shaped transducer on a plane using the Rayleigh integral. The integral
% is calculated at a distance between the transducer’s edge and the 
% focus point. The resulting field is then propagated to the desired 
% output plane using the angular spectrum approach (ASA). The output 
% plane can either coincide with the center of the bowl (z = 0) or be
% moved further to adjust the size of the simulation volume. The ASA 
% propagation technique is adapted from the k-Wave toolbox. 
%
% INPUTS:
% 
% c0: Sound speed in the medium [m/s]
% rho0: Medium density [kg/m³]
% a: Aperture radius of the bowl-shaped transducer [m]
% r: Radius of curvature of the bowl [m]
% f: Source frequency [Hz]
% dx: Grid spacing in both input and output planes [m]
% v0: oscillation velocity on the bowl surface [m/s]
% source_pos: Position of the hologram plane. A value of 0 coincides with 
% the center of the bowl’s surface, while positive values shift along 
% the beam propagation direction [m].

% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%     'Padding'      - Grid padding used to increase the accuracy of
%                      the projection. Can be modified if the accuracy is 
%                      not enough (default = 0).
% define defaults
padding      = 0;

% replace with user defined values if provided
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Padding'                
                % assign input
                padding = round(varargin{input_index + 1});                          
            otherwise
                error('Unknown optional input.');
        end
    end
end
                      
pml_layer = 40; % usual pml layer thickness is 20 point from both sides
min_size_x = pml_layer + round(1.4 * a * 2/dx); %points in plane size, [m]/grid step

%find an optimal size of grid to make an odd amount of points and low prime
%factor
Nx = find_odd_factor(min_size_x, min_size_x+40);

%find a suitable plane for Rayleigh integral calculation - several
%wavelengths from the bowl edge
bowl_thickness = r - sqrt(r^2 - a^2); %thickness of the bowl
wavelength = c0/f; %wavelength
dist_for_calc = bowl_thickness + 3 * wavelength; % position of 
%the hologram plane for Rayleigh integral calculation

%Now set coordinates of points at the computational plane
[x_hologram, y_hologram] = meshgrid((-(Nx+1)/2 : (Nx+1)/2)*dx);
z_hologram = dist_for_calc * ones(size(x_hologram));

%Set coordinates of point at the surface of the bowl
x_vector_source = linspace(-a,a, fix(2*a/(0.3*c0/f)*2));
y_vector_source = linspace(-a,a, fix(2*a/(0.3*c0/f)*2));

[y_source, x_source] = meshgrid(x_vector_source,y_vector_source);
z_source = r - sqrt(r^2 - x_source.^2 - y_source.^2);

%Set initial velocity on the surface
VelocityInitial = zeros(size(x_source));
VelocityInitial(x_source.^2 + y_source.^2 <= a^2) = v0;

loop_start_time = clock;

%Compute Raileigh integral on a simulation plane
prHolo = CalcRay( x_source, y_source, z_source, x_hologram, y_hologram, z_hologram, VelocityInitial, rho0, c0, f, -1, 'CalcSphericalSourceInitialV', r);
%Transfer the hologram to the desired plane at source_pos by ASA
prHolo_source = angularSpectrumCW_edited(prHolo, dx, source_pos - dist_for_calc, f, c0, 'GridExpansion', padding);
%Compute velocity hologram to use in k_wave
vHolo = conj(holoP_2_holoV(prHolo_source, dx, dx, f, c0, rho0));

disp(['Hologram calculation completed in ' scaleTime(etime(clock, loop_start_time))]);

%TO SET AS A SOURCE IN K-WAVE DO THE FOLLOWING:
% ampl = abs(vHolo;
% phase = angle(vHolo);% 
% source.p = (sin(2*pi*TRANS_FREQ*(0:(kgrid.Nt - 1)) * dt).*cos(phase(:)) - cos(2*pi*TRANS_FREQ*(0:(kgrid.Nt - 1)) * dt).*sin(phase(:))).*ampl(:);
    
    
%Check the accuracy - axial profiles must coincide well outside the near-field   
    
%calculate teh axial distribution 
[p_axial, ~] = focusedBowlONeil(r, ...
          2*a, v0, f, c0, rho0, ...
          0:dx:1.5*r, 0);

pr_ax_asa = angularSpectrumCW_edited(prHolo_source, dx, 0:dx:1.5*r, f, c0, 'GridExpansion', 0);
[~, ind] = maxND(pr_ax_asa);

figure;
plot(1e3*(0:dx:1.5*r), p_axial, 'black', 'LineWidth', 2);
hold on;
plot( 1e3*(0:dx:1.5*r), abs(squeeze(pr_ax_asa(ind(1),ind(2),:))), 'r--', 'LineWidth', 2)
hold on
legend('Analytical solution', 'ASA expansion of hologram');
title('Make sure plots coincide well in far-field!');
xlabel('Z [mm]');
ylabel('P [Pa]');

disp(['Accuracy check completed in ' scaleTime(etime(clock, loop_start_time))]);

end



function pressure = angularSpectrumCW_edited(input_plane, dx, z_pos, f0, medium, varargin)
% ANGULARSPECTRUMCW Project CW input plane using the angular spectrum method.
%
% DESCRIPTION:
%     EDITED - commented the last section % trim grid expansion %

%     angularSpectrumCW projects an input plane of single-frequency
%     continuous wave data (given as a 2D matrix of complex pressure
%     values) to the parallel plane or planes specified by z_pos using the
%     angular spectrum method. The implementation follows the spectral
%     propagator with angular restriction described in reference [1].
%
%     For linear projections in a lossless medium, just the sound speed can
%     be specified. For projections in a lossy medium, the parameters are
%     given as fields to the input structure medium. 
%
%     To compute the pressure field over an isotropic domain with Nz grid
%     points (assuming the source plane is aligned with z_ind = 1), use the
%     syntax: 
%
%           pressure = angularSpectrumCW(input_plane, dx, (0:(Nz - 1)) * dx, f0, c0)
%
% USAGE:
%     pressure = angularSpectrumCW(input_plane, dx, z_pos, f0, c0)
%     pressure = angularSpectrumCW(input_plane, dx, z_pos, f0, c0, ...)
%     pressure = angularSpectrumCW(input_plane, dx, z_pos, f0, medium)
%     pressure = angularSpectrumCW(input_plane, dx, z_pos, f0, medium, ...)
%
% INPUTS:
%     input_plane          - 2D matrix of complex pressure values over a
%                            plane [Pa].
%     dx                   - Spatial step between grid points in the input
%                            plane [m].
%     z_pos                - Vector specifying the relative z-position of
%                            the planes to which the data is projected [m].
%     f0                   - Source frequency [Hz].
%
%     c0                   - Medium sound speed [m/s].
%             OR
%     medium.sound_speed   - Medium sound speed [m/s].
%     medium.alpha_power   - Power law absorption exponent.
%     medium.alpha_coeff   - Power law absorption coefficient
%                            [dB/(MHz^y cm)].
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'AngularRestriction' - Boolean controlling whether angular
%                            restriction is used as described in [1]
%                            (default = true).
%     'DataCast'           - String input of the data type that variables
%                            are cast to before computation. For example,
%                            setting to 'single' will speed up the
%                            computation time (due to the improved
%                            efficiency of fft2 and ifft2 for this data
%                            type). This variable is also useful for
%                            utilising GPU parallelisation the Parallel
%                            Computing Toolbox by setting 'DataCast' to
%                            'gpuArray-single' (default = 'off').
%     'DataRecast'         - Boolean controlling whether the output data
%                            is cast back to double precision. If set to
%                            false, sensor_data will be returned in the
%                            data format set using the 'DataCast' option
%                            (default = false).
%     'FFTLength'          - Length of the FFT used to compute the angular
%                            spectrum (default = 1 + the next power of two
%                            larger than the grid size).
%     'GridExpansion'      - Grid padding used to increase the accuracy of
%                            the projection. The grid expansion is removed
%                            before returning the calculated pressure to
%                            the user (default = 0).
%     'Reverse'            - Boolean controlling whether the projection is
%                            in the forward (false) or backward (true)
%                            direction (default = false).
%
% OUTPUTS:
%     pressure             - 3D matrix of complex pressure values across
%                            the 2D planes specified by z_pos, indexed as
%                            (x_ind, y_ind, plane_index) [Pa].
% 
% ABOUT:
%     author               - Bradley Treeby
%     edited               - Alisa Krokhmal 27th September 2024
%     date                 - 27th February 2018
%     last update          - 19th February 2019
%
% REFERENCES:
%     [1] Zeng, X., & McGough, R. J. (2008). Evaluation of the angular
%     spectrum approach for simulations of near-field pressures. The
%     Journal of the Acoustical Society of America, 123(1), 68-76.
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby
%
% See also angularSpectrum

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

% start timer
start_time = clock;

% =========================================================================
% INPUT CHECKING
% =========================================================================

% define defaults
angular_restriction = true;
grid_expansion      = 0;
fft_length          = 'auto';
data_cast           = 'off';
data_recast         = false;
reverse_proj        = false;
absorbing           = false;
loops_for_time_est  = 5;

% replace with user defined values if provided
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'AngularRestriction'
                
                % assign input
                angular_restriction = logical(varargin{input_index + 1});
                
            case 'DataCast'
                
                % assign input
                data_cast = varargin{input_index + 1};
                
                % check list of valid inputs
                if ~ischar(data_cast)
                    error('Optional input ''DataCast'' must be a string.');
                elseif ~(strcmp(data_cast, 'off') || strcmp(data_cast, 'double') ...
                        || strcmp(data_cast, 'single') || strcmp(data_cast, 'gpuArray-single') ... 
                        || strcmp(data_cast, 'gpuArray-double'))
                    error('Invalid input for ''DataCast''.');
                end
                
                % replace double with off
                if strcmp(data_cast, 'double')
                    data_cast = 'off';
                end
                
                % create empty string to hold extra cast variable for use
                % with the parallel computing toolbox
                data_cast_prepend = '';
                
                % replace PCT options with gpuArray
                if strcmp(data_cast, 'gpuArray-single')
                    data_cast = 'gpuArray';
                    data_cast_prepend = 'single';
                elseif strcmp(data_cast, 'gpuArray-double')
                    data_cast = 'gpuArray';
                end
                
                if strcmp(data_cast, 'gpuArray')
                    
                    % check the PCT is installed and the version is 2012a or
                    % later (verLessThan only works in release 7.4 or later)
                    v = ver;
                    if verLessThan('matlab', '7.14') || ~ismember('Parallel Computing Toolbox', {v.Name})
                        error('The Parallel Computing Toolbox for MATLAB 2012a or later is required for ''DataCast'' set to ''gpuArray-single'' or ''gpuArray-double''.');
                    end

                    % cleanup unused variables
                    clear v;
                    
                end
                
            case 'DataRecast'
                
                % assign input
                data_recast = logical(varargin{input_index + 1});
                
            case 'FFTLength'
                
                % assign input
                fft_length = round(varargin{input_index + 1});
                
            case 'GridExpansion'
                
                % assign input
                grid_expansion = round(varargin{input_index + 1});
                
            case 'Reverse'
                
                % assign input
                reverse_proj = logical(varargin{input_index + 1});
                
            otherwise
                error('Unknown optional input.');
        end
    end
end

% check for structured medium input
if isstruct(medium)
    
    % force the sound speed to be defined
    if ~isfield(medium, 'sound_speed')
        error('medium.sound_speed must be defined when specifying medium properties using a structure.');
    end
    
    % assign the sound speed
    c0 = medium.sound_speed;
    
    % assign the absorption
    if isfield(medium, 'alpha_coeff') || isfield(medium, 'alpha_power')
        
        % enforce both absorption parameters
        if ~(isfield(medium, 'alpha_coeff') && isfield(medium, 'alpha_power'))
            error('Both medium.alpha_coeff and medium.alpha_power must be defined for an absorbing medium');
        end
        
        % convert attenuation to Np/m
        alpha_Np = db2neper(medium.alpha_coeff, medium.alpha_power) * (2 * pi * f0)^medium.alpha_power;
        
        % check for zero absorption and assign flag
        if alpha_Np ~= 0
            absorbing = true;
        end
        
    end
    
else
    
    % assign the sound speed
    c0 = medium;
    
end

% check for maximum supported frequency
if dx > (c0/(2 * f0))
    error(['Input frequency is higher than maximum supported frequency of ' scaleSI(c0 / (2 * dx)) 'Hz.' ]);
end

% =========================================================================
% PRE-PROCESSING
% =========================================================================

% get grid size
[Nx, Ny] = size(input_plane);
Nz = length(z_pos);

% get scale factor for grid size
[~, scale, prefix] = scaleSI(min(Nx * dx, Ny * dx));

% update command line status
% disp('Running CW angular spectrum projection...');
% disp(['  start time: ' datestr(start_time)]);
% disp(['  input plane size: ' num2str(Nx) ' by ' num2str(Ny) ' grid points (' num2str(scale * Nx * dx) ' by ' num2str(scale * Ny * dx) prefix 'm)']);
% disp(['  grid expansion: ' num2str(grid_expansion) ' grid points']);

% apply phase conjugation if stepping backwards
if reverse_proj
    input_plane = conj(input_plane);
end

% expand input
if grid_expansion > 0
    input_plane = expandMatrix(input_plane, grid_expansion, 0);
    [Nx, Ny] = size(input_plane);
end

% get FFT size
if ischar(fft_length) && strcmp(fft_length, 'auto')
    fft_length = 2.^(nextpow2(max(size(input_plane))) + 1);
end

% update command line status
% disp(['  FFT size: ' num2str(fft_length) ' points']);
% disp(['  maximum supported frequency: ' scaleSI( c0 / (2 * dx) ) 'Hz']);

% create wavenumber vector
N = fft_length;
if rem(N, 2) == 0
    k_vec = ((-N/2):(N/2-1)) .* 2 * pi ./ (N * dx);
else
    k_vec = (-(N-1)/2:(N-1)/2) .* 2 * pi ./ (N * dx);
end

% force middle value to be zero in case 1/Nx is a recurring
% number and the series doesn't give exactly zero
k_vec(floor(N/2) + 1) = 0;

% shift wavenumbers to be in the correct order for FFTW
k_vec = ifftshift(k_vec);

% compute wavenumber at driving frequency
k = 2 * pi * f0 / c0;

% create wavenumber grids
[kx, ky] = meshgrid(k_vec, k_vec);
kz = sqrt(k.^2 - complex(kx.^2 + ky.^2));
 
% precompute term for angular restriction
sqrt_kx2_ky2 = sqrt(kx.^2 + ky.^2);

% preallocate output and copy source plane
pressure = zeros(Nx, Ny, Nz);
pressure(:, :, 1) = input_plane;

% compute forward Fourier transform of input plane
input_plane_fft = fft2(input_plane, fft_length, fft_length);

% =========================================================================
% DATA CASTING
% =========================================================================

% cast variables to the output type
if ~strcmp(data_cast, 'off')
    
    % update command line status
    disp(['  casting variables to ' data_cast ' type...']);    
    
    % list of variables to cast
    cast_variables = {'kz', 'z_pos', 'input_plane_fft', 'pressure'};
    
    % additional variables used if absorbing
    if absorbing
        cast_variables = [cast_variables, {'alpha_Np', 'k'}];
    end
        
    % additional variables used for angular restriction
    if angular_restriction
        cast_variables = [cast_variables, {'sqrt_kx2_ky2', 'fft_length', 'dx'}];
    end
    
    % loop through, and change data type
    for cast_index = 1:length(cast_variables)
        eval([cast_variables{cast_index} ' = ' data_cast '(' data_cast_prepend '(' cast_variables{cast_index} '));']);
    end
    
end

% =========================================================================
% Z-LOOP
% =========================================================================

% update command line status
loop_start_time = clock;
% disp(['  precomputation completed in ' scaleTime(etime(loop_start_time, start_time))]);
% disp('  starting time loop...');

% loop over z-positions
for z_index = 1:Nz

    % get current z value
    z = z_pos(z_index);

    % if set to zero, just store the input plane
    if z == 0
        
        % store input data
        pressure(:, :, z_index) = input_plane;
        
    else
    
        % compute spectral propagator (Eq. 6)
        H = conj(exp(1i .* z .* kz));

        % account for attenuation (Eq. 11)
        if absorbing
            H = H .* exp(-alpha_Np .* z .* k ./ kz);
        end

        % apply angular restriction
        if angular_restriction

            % size of computational domain [m]
            D = (fft_length - 1) * dx;

            % compute angular threshold (Eq. 10)
            kc = k * sqrt(0.5 * D.^2 ./ (0.5 * D.^2 + z.^2));

            % apply threshold to propagator
            H(sqrt_kx2_ky2 > kc) = 0;

        end

        % compute projected field and store
        pressure_step = ifft2(input_plane_fft .* H, fft_length, fft_length);
        pressure(:, :, z_index) = pressure_step(1:Nx, 1:Ny);
        
    end
    
    % update command line status
%     if z_index == loops_for_time_est
%         disp(['  estimated simulation time ' scaleTime(etime(clock, loop_start_time) * Nz / z_index) '...']);
%     end
    
end

% update command line status
% disp(['  simulation completed in ' scaleTime(etime(clock, loop_start_time))]);

% =========================================================================
% POST PROCESSING
% =========================================================================

% trim grid expansion
% if grid_expansion > 0
%     pressure = pressure(1 + grid_expansion:end - grid_expansion, 1 + grid_expansion:end - grid_expansion, :);
% end

% reverse grid and conjugate the phase
if reverse_proj
    pressure = flip(conj(pressure), 3);
end

% cast output back to double precision
if data_recast
    if strncmp(data_cast, 'gpuArray', 8)
        pressure = double(gather(pressure));
    else
        pressure = double(pressure);
    end
end

% update command line status
% disp(['  total computation time ' scaleTime(etime(clock, start_time))]);
end

function min_factor = find_odd_factor(min_number, max_number)
% Return minumum odd number from the range with whe minimum prime factor
%
% DESCRIPTION:
%     find_odd_factor loops through the given range of numbers and finds the
%     odd number with the smallest maximum prime factors. This allows suitable
%     grid sizes to be selected to maximise the speed of the FFT. Odd number is
%     necessary for the field calculation of the axial symmetric source to compare 
%     axial profiles. The prime factor is printed to the command line.
%    
%
% INPUTS:
%     min_number    - integer specifying the lower bound of values to test
%     max_number    - integer specifying the upper bound of values to test
%
% OUTPUT:
%     min_factor    - integer odd number from the given range with the minimum prime factor
%
%
% ABOUT:
%     author        - Alisa Krokhmal
%     last update   - 26th August 2024

% extract factors
nums = min_number:max_number;
odd_nums = nums(mod(nums,2)==1);

facs = zeros(size(odd_nums));
fac_max = facs;
min_odd = min(odd_nums);

for index = 1:length(odd_nums)
    facs(index) = length(factor(odd_nums(index)));
    fac_max(index) = max(factor(odd_nums(index)));
end

% compute best factors in range
disp('Prime factor is');
ind = min(fac_max);
disp(num2str(ind));

min_factor =  odd_nums(find(fac_max == min(fac_max), 1));
end

function [ field_out ] = CalcRay( x_source, y_source, z_source, x_hologram, y_hologram, z_hologram, field_in, rho0, c0, f0, exp_sign, type, Focal_distance)

%DESCRIPTION:
% CalcRay computes the Rayleigh integral from the set type of the source on
% a plane or on a sphere
%
%author: Pavel Rosnitskiy

  if nargin < 13
    Focal_distance =   Inf;
  end
  
percent_step = 1;
i_percent_step = 1;


k = 2*pi*f0/c0;
i_compl = exp_sign*1i;



if strcmp(type,'BackPlaneSourceP2V')||strcmp(type,'BackSphericalSourceP2V')
    dx = x_hologram(2,2) - x_hologram(1,1);
    dy = y_hologram(2,2) - y_hologram(1,1);
    
    field_row = zeros(1,numel(x_source));
    for i1 = 1:numel(x_source)
        
        
        percent = i1/numel(x_source)*100;
        if  (percent > percent_step*i_percent_step)
            disp(['Calculation of the Rayleigh integral...' num2str(fix(percent)) '% done']);
            i_percent_step = i_percent_step + 1;
        end
        R = sqrt((x_hologram-x_source(i1)).^2 + (y_hologram-y_source(i1)).^2 + (z_hologram-z_source(i1)).^2);
        
        if strcmp(type,'BackPlaneSourceP2V')
            %             K = 1/2/pi/i_compl/k/rho0/c0*(-1*(i_compl*k./R + R.^-2) +  ((z_hologram-z_source(i1))./R).^2.*(3*i_compl*k./R + 3./R.^2 - k^2)).*exp(-i_compl*k*R)./R;
            n1n2  = -1;
            m12n1 = (z_hologram-z_source(i1))./R;
            m21n2 = (z_hologram-z_source(i1))./R;
            
        elseif strcmp(type,'BackSphericalSourceP2V')
            
            n1n2  = z_source(i1)/Focal_distance - 1;
            m12n1 = 1./R.*((Focal_distance^-1*(x_source(i1)^2 + y_source(i1)^2 + z_source(i1)^2 - x_source(i1)*x_hologram - y_source(i1)*y_hologram - z_source(i1)*z_hologram)) + z_hologram-z_source(i1));
            m21n2 = (z_hologram-z_source(i1))./R;
            
        end
        K = 1/2/pi/i_compl/k/rho0/c0*(n1n2*(i_compl*k./R + R.^-2) +  m12n1.*m21n2.*(3*i_compl*k./R + 3./R.^2 - k^2)).*exp(-i_compl*k*R)./R;
        field_row(i1)    = sum(sum(field_in.*K))*dx*dy;
        
    end
    field_out = reshape(field_row, size(x_source));
else
    dx = x_source(2,2) - x_source(1,1);
    dy = y_source(2,2) - y_source(1,1);
    
    field_row = zeros(1,numel(x_hologram));
    for i2 = 1:numel(x_hologram)
        
        
        percent = i2/numel(x_hologram)*100;
        if  (percent > percent_step*i_percent_step)
            disp(['Calculation of the Rayleigh integral...' num2str(fix(percent)) '% done']);
            i_percent_step = i_percent_step + 1;
        end
        R = sqrt((x_hologram(i2)-x_source).^2 + (y_hologram(i2)-y_source).^2 + (z_hologram(i2)-z_source).^2);
        if strcmp(type,'CalcSphericalSourceInitialV')
            K = -i_compl*k*rho0*c0/2/pi.*exp(i_compl*k*R)./R./(z_source/Focal_distance - 1);
        elseif strcmp(type,'CalcPlaneSourceInitialV')
            K = -i_compl*k*rho0*c0/2/pi.*exp(i_compl*k*R)./R;
        elseif strcmp(type,'CalcPlaneSourceInitialP')
            K = 1/2/pi*(z_hologram(i2)-z_source)./R.*(-i_compl*k./R + R.^-2).*exp(i_compl*k*R);
        end
        field_row(i2)    = sum(sum(field_in.*K))*dx*dy;
        
    end
    field_out = reshape(field_row, size(x_hologram));
end

end

function [V] = holoP_2_holoV(P, dx, dy, f0, c0, rho0)
% Conversion of the pressure hologram to velocity hologram to be suitable
% to set as a source in k-Wave
%Author: Pavel Rosnitskiy


k = 2*pi*f0/c0;

if (size(P,1) == size(P,2))
    Nx = size(P,1);
    Lx = (Nx-1)*dx;
    
    Ny = size(P,1);
    Ly = (Ny-1)*dy;
else
    error('Pressure matrix is not square!')
end

kx = -(pi/dx):(2*pi/Lx):(pi/dx);
ky = -(pi/dy):(2*pi/Ly):(pi/dy);
[Kx, Ky] = meshgrid(kx, ky);

Sp = fft2(P);
Sp = fftshift(Sp);

Sv = 1/rho0/c0*sqrt(1-(Kx.^2 + Ky.^2)/k^2).*Sp;

Sv = ifftshift(Sv);

V = ifft2(Sv);

end
