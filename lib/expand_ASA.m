function p_all = expand_ASA(p, p_max, dt, dx, z_sensor, z_lim,  TRANSDUCER_FREQ, WATER_TEMPERATURE)

    FILTER_FREQ = [0.1e6 3e6 0.1e6]; % Hz, [start, end, taper width]
    PROJ_DISTANCE_BACK = z_sensor; %%%%%%%%%ЗАДАЙ ЗДЕСЬ 30 для шага 0.5 или 12.5 для шага 0.125!
    PROJ_DISTANCE_FORWARD = z_lim;

    p_max = permute(p_max, [2, 3, 1]);
    N_shape = size(p);
    Nxy = (sqrt(N_shape(1)));

    % calculate Nz based on the point spacing dx, and the distance to the
    % transducer 
    c0 = waterSoundSpeed(WATER_TEMPERATURE);
    Nz_back = round(PROJ_DISTANCE_BACK / dx); 
    Nz_forward = round(PROJ_DISTANCE_FORWARD / dx); 

    [mag, phase] = extractAmpPhase(p, 1/dt, TRANSDUCER_FREQ, 'Dim', 2, 'FFTPadding', 3);
    size_p = size(p);
    p_plane= reshape(mag.*exp(1i*phase), sqrt(size_p(1)), sqrt(size_p(1)));

    pressure_forward = abs(angularSpectrumCW(p_plane, dx, 0:dx:PROJ_DISTANCE_FORWARD, TRANSDUCER_FREQ, c0, 'Reverse', 0, 'GridExpansion', 50));
    z_total = (0:dx:PROJ_DISTANCE_FORWARD+PROJ_DISTANCE_BACK);


    p_all = cat(3, p_max(:,:, 21:21+round(PROJ_DISTANCE_BACK/dx)-1), pressure_forward(:, :, :));
end
% [~, ind] = maxND(pressure_forward);

% figure;
% imagesc(z_total*1e3, (-Nxy/2:Nxy/2)*dx*1e3, squeeze(p_all(ind(1), :, :)));
% figure;
% plot(z_total*1e3, squeeze(p_all(ind(1), ind(2), :))/1e3);






 
 

