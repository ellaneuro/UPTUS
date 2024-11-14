function UPTUS_sim_3D(t1_filename, ct_filename, output_dir, focus_coords_in, bowl_coords_in, focus_depth, ppw, run_simulation_ac, run_thermal_sim, pulse_dur, pulse_rep_int, pulse_train_dur, pressure, phase, transducer)

%% = change with each run/entry point - add to config

% transducer info DPX-500
% element amplitude 20000
% element outer diameter 0.0322
% element inner diameter 0.0054
% elements type 1
% driving frequency 50000
% target pressure 5.4772e+05
% focal distance 0.1000
% focal distance max 0.1480
% focal distance min 0.0360
% focal distance surface offset 0.0051
% target pressure default 5.4772e+05
% target pressure max 9.4868e+05

run_simulation_ac = true;
run_thermal_sim = false;

t1_filename = '/Users/ellasonnelid/UPTUS/input/DICOM_sT1W_3D_TFE_20241022144257_501.nii'; %%
ct_filename = '/Users/ellasonnelid/UPTUS/output/sub_001_pct.nii'; %%
output_dir = '/Users/ellasonnelid/UPTUS/output'; %% använd ny för sim? output_sim?
focus_coords_in = [121, 143, 116]; %% add 1 to all for Matlab
bowl_coords_in = [105, 210, 146]; %% add 1 to all for Matlab
focus_depth = 75; %%

% entry_point_ras = [42.51, 83.62, 82.41] %% byt ut mot voxel space direkt? 
% target_point_ras = [-7.76, 29.75, 44.20]

ppw = 6; % changed to 2 instead of 3 in attempt to make the simulations less heavy 
pulse_dur = 20e-3
pulse_rep_int = 200e-3
pulse_train_dur = 120e-3
pressure = 70000 % p = (I⋅ρ⋅c)^(1/2), I = 30 w/cm^2, p = density of brain tissue, c = speed of sound in brain %% dubbel kolla! 
phase = [0, 104, 209, 145] %%

% medium parameter
% https://pmc.ncbi.nlm.nih.gov/articles/PMC8891811/
% ITRUSST benchmarks - https://github.com/Donders-Institute/PRESTUS/blob/master/configs/default_config.yaml
c_min               = 1500;     % sound speed [m/s]
c_max               = 3100;     % max. speed of sound in skull (F. A. Duck, 2013.) [m/s]
rho_min             = 1000;     % density [kg/m^3]
rho_max             = 1900;     % max. skull density [kg/m3]
alpha_power         = 1.43;     % Robertson et al., PMB 2017 usually between 1 and 3? from Treeby paper
alpha_coeff_water   = 0;        % [dB/(MHz^y cm)] close to 0 (Mueller et al., 2017), see also 0.05 Fomenko et al., 2020?
alpha_coeff_min     = 4;        
alpha_coeff_max     = 8.7; 

% transducer DPX-500
transducer = struct();
transducer.roc = 0.150; 
transducer.diameters = [0.0054 0.0322; 0.0322 0.0454; 0.0454 0.0555; 0.0555 0.0640];
transducer.diameters = reshape(transducer.diameters', 2, []);
transducer.frequency = 500e3;

record_periods  = 3;
cfl             = 0.3;

hu_min 	= 300;
hu_max 	= 2000; 

input_ct = niftiread(ct_filename);
input_t1 = niftiread(t1_filename);
header = niftiinfo(ct_filename);

dx = c_min / (ppw * transducer.frequency);
scale_factor = 0.1 % round(header.PixelDimensions/(dx*1e3),2); % took this out to see if it make the size smaller 
ct_img = imresize3(input_ct, 'cubic', 'Scale', scale_factor);
t1_img = imresize3(input_t1, 'cubic', 'Scale', scale_factor);

focus_coords = round(focus_coords_in.*scale_factor);
bowl_coords = round(bowl_coords_in.*scale_factor);

ct_max = max(ct_img(:));
if ct_max < hu_max
    hu_max = ct_max;
end
clear ct_max;

skull_model = ct_img;
skull_model(skull_model < hu_min) = 0;
skull_model(skull_model > hu_max) = hu_max;

midpoint = round((bowl_coords + focus_coords)/2);
padx = 180;
tmp_model = zeros(size(skull_model,1)+padx*2, size(skull_model,2)+padx*2, size(skull_model,3)+padx*2);
tmp_model(padx+1:size(skull_model,1)+padx, ...
    padx+1:size(skull_model,2)+padx, ...
    padx+1:size(skull_model,3)+padx) = skull_model;
tmp_midpoint = midpoint+padx;

shift_idx = [tmp_midpoint(1)-padx+1,tmp_midpoint(1)+padx;
    tmp_midpoint(2)-padx+1,tmp_midpoint(2)+padx; ...
    tmp_midpoint(3)-padx+1,tmp_midpoint(3)+padx];
model = tmp_model(shift_idx(1,1):shift_idx(1,2), ...
    shift_idx(2,1):shift_idx(2,2), ...
    shift_idx(3,1):shift_idx(3,2));

shift_x = padx-midpoint(1);
shift_y = padx-midpoint(2);
shift_z = padx-midpoint(3);

bowl_coords     = bowl_coords + [shift_x, shift_y, shift_z];
focus_coords    = focus_coords + [shift_x, shift_y, shift_z];

new_t1 = zeros(size(model));
idx1 = zeros(3,2);
for ii = 1:3
    if shift_idx(ii,1)-padx <= 0
        idx1(ii,1) = 1;
    elseif shift_idx(ii,1)-padx > 0
        idx1(ii,1) = shift_idx(ii,1)-padx;
    end
    if shift_idx(ii,2)-padx <= size(ct_img,ii)
        idx1(ii,2) = shift_idx(ii,2)-padx;
    elseif shift_idx(ii,2)-padx > size(ct_img,ii)
        idx1(ii,2) = size(ct_img,ii);
    end
end
idx2 = zeros(3,2);
idx2(:,1) = [shift_x, shift_y, shift_z] + [1,1,1];
idx2(:,2) = size(ct_img) + [shift_x, shift_y, shift_z];
for ii = 1:3
    if idx2(ii,1) <= 0
        idx2(ii,1) = 1;
    end
    if idx2(ii,2) > size(model,ii)
        idx2(ii,2) = size(model,ii);
    end
end
new_t1(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2)) = ...
    t1_img(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2), idx1(3,1):idx1(3,2));
t1_img = new_t1;
clear new_t1 tmp_model;

medium.density = rho_min + (rho_max - rho_min) * ...
    (model - 0) / (hu_max - 0);
medium.sound_speed = c_min + (c_max - c_min) * ...
    (medium.density - rho_min) / (rho_max - rho_min);
medium.alpha_coeff = alpha_coeff_min + (alpha_coeff_max - alpha_coeff_min) * ...
    (1 - (model - hu_min) / (hu_max - hu_min)).^0.5;

medium.density(model == 0) = rho_min;
medium.sound_speed(model == 0) = c_min;
medium.alpha_coeff(model == 0) = alpha_coeff_water;

medium.alpha_power = alpha_power;

[Nx, Ny, Nz] = size(model);
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

ppp = round(ppw / cfl);

dt = 1 / (ppp * transducer.frequency);
dt_stability_limit = checkStability(kgrid, medium);
if dt_stability_limit ~= Inf
    dt = dt_stability_limit;
end

t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2) / c_min;

Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

cfl = c_min * dt / dx;
ppp = round(ppw / cfl);
disp(['PPW = ' num2str(ppw)]);
disp(['CFL = ' num2str(cfl)]);
disp(['PPP = ' num2str(ppp)]);

bowl_pos = [kgrid.x_vec(bowl_coords(1)), ...
    kgrid.y_vec(bowl_coords(2)), ...
    kgrid.z_vec(bowl_coords(3))];
focus_pos = [kgrid.x_vec(focus_coords(1)), ...
    kgrid.y_vec(focus_coords(2)), ...
    kgrid.z_vec(focus_coords(3))];

karray = kWaveArray('SinglePrecision', true);
karray.addAnnularArray(bowl_pos, transducer.roc, transducer.diameters, focus_pos)

source.p_mask = karray.getArrayBinaryMask(kgrid);
if ~run_simulation_ac && ~run_thermal_sim
    model_mask = model;
    model_mask(model_mask>0)=1;
    tmp = source.p_mask+model_mask;
    volumeViewer(tmp);
    clear model_mask tmp;
    return
end


if run_simulation_ac

    source_amp = repmat(pressure, [1, 4]);
    source_phase = deg2rad(phase);
    freq = transducer.frequency;

    % time varying source
    source_sig = createCWSignals(kgrid.t_array, freq, source_amp, source_phase);
    source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

    % sensor
    sensor.mask = ones(Nx, Ny, Nz);
    sensor.record = {'p'};
    sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

    input_args = {'PMLSize', 10, 'PMLInside', true, 'PlotPML', true, ... 
    'DisplayMask', source.p_mask};

    if ~exist(fullfile(output_dir), 'dir')
        mkdir(output_dir)
    end

    input_filename = fullfile(output_dir, 'kwave_input.h5');
    output_filename = fullfile(output_dir, 'kwave_output.h5');

    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

    % get sensor data - amplitude
    p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

    p = reshape(p, Nx, Ny, Nz);

    [max_pressure, idx] = max(p(:));  % [Pa]
    [mx, my, mz] = ind2sub(size(p), idx);

    % max Isppa
    Isppa = max_pressure^2 / (2 * rho_min * c_min);
    Isppa = Isppa * 1e-4;

    % mechanical index
    MI = max_pressure * 1e-6 / sqrt(freq * 1e-6);

    tmp_focal_vol = int16(p > 0.71 * max(p(:)));
    cc = bwconncomp(tmp_focal_vol);
    focal_vol = length(cc.PixelIdxList{1}) * (dx * 1e3)^3;
    clear tmp_focal_vol cc;

    % pressure + Isppa at target 
    p_focus = p(focus_coords(1), focus_coords(2), focus_coords(3));
    isppa_focus = p_focus^2 / (2 * rho_min * c_min) * 1e-4; 

    if norm(focus_coords-[mx,my,mz])*dx*1e3 > 5
        warning(['Maximum pressure point is more than 5 mm away from the intended focus. ' ...
            'Attempting to check whether max pressure point is at the skull interface.'])

        u = ([mx,my,mz]-focus_coords)/norm([mx,my,mz]-focus_coords);
        tmp_max = round([mx,my,mz] - 2*u);

        if model(tmp_max(1), tmp_max(2), tmp_max(3)) > 1
            warning(['It is likely that the maximum pressure is at the skull interface. ' ...
                'Attempting to adjust search for maximum pressure along the trajectory... ' ...
                'I cannot guarantee this will be correct... ' ...
                'Saving .mat file, please check output!'])
            clear source sensor_data;
            save(fullfile(output_dir, 'acoustic_sim_output.mat'));

            u = ([mx,my,mz]-focus_coords)/norm([mx,my,mz]-focus_coords);
            search_rad = 10;
            search_start = [mx,my,mz] - search_rad*u;
            search_start = round(search_start);
            tmp=p(search_start(1):end,1:search_start(2),1:search_start(3));

            [max_pressure, idx] = max(tmp(:));
            [mx, my, mz] = ind2sub(size(tmp), idx);

            mx = search_start(1) + mx;
            my = search_start(2) - my;
            mz = search_start(3) - mz;

            Isppa = max_pressure^2 / (2 * rho_min * c_min);
            Isppa = Isppa * 1e-4;

            MI = max_pressure * 1e-6 / sqrt(freq * 1e-6);

            tmp_focal_vol = int16(p>0.71*max_pressure);
            cc = bwconncomp(tmp_focal_vol);
            focal_vol = length(cc.PixelIdxList{1})*(dx*1e3)^3;
            clear tmp_focal_vol cc;

        elseif model(tmp_max(1), tmp_max(2), tmp_max(3)) == 0
            warning(['Maximum pressure is probably not at the skull interface, ' ...
                'but it is still more than 5 mm away from your planned focus... ' ...
                'Please check output!'])
        end
    end
    
    figure;
    ax1 = axes; imagesc(ax1, imrotate(squeeze(t1_img(mx,:,:)),90), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, imrotate(squeeze(p(mx,:,:))*1e-6,90));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'turbo')
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
    xlabel(cb2, '[MPa]');
    title(ax1,'Acoustic Pressure Amplitude')
    saveas(gcf, fullfile(output_dir, 'pressure_sag.jpg'));

    figure;
    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_img(mx,:,:)),90), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, imrotate(squeeze(p(mx,:,:)>(0.71*max_pressure))*1e-6,90));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'turbo')
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
    xlabel(cb2, '[MPa]');
    title(ax1,'50% Acoustic Pressure Amplitude')
    saveas(gcf, fullfile(output_dir, 'pressure_sag_50%.jpg'));

    figure;
    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_img(:,my,:)),90), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, imrotate(squeeze(p(:,my,:))*1e-6,90));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'turbo')
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
    xlabel(cb2, '[MPa]');
    title(ax1,'Acoustic Pressure Amplitude')
    saveas(gcf, fullfile(output_dir, 'pressure_cor.jpg'));

    figure;
    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_img(:,:,mz)),90), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, imrotate(squeeze(p(:,:,mz))*1e-6,90));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'turbo')
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
    xlabel(cb2, '[MPa]');
    title(ax1,'Acoustic Pressure Amplitude')
    saveas(gcf, fullfile(output_dir, 'pressure_ax.jpg'));

    p_out = zeros(size(ct_img),'single'); % changed to single in attemp to make simulations less heavy
    p_out(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2), idx1(3,1):idx1(3,2)) = ...
        p(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2));
    p_out = imresize3(p_out, 'cubic', 'Scale', 1./scale_factor);
    header.ImageSize = size(p_out);
    header.Filename=[]; header.Filemoddate=[]; header.Filesize=[]; header.raw=[];
    header.Datatype='double'; header.BitsPerPixel=32;
    niftiwrite(p_out, fullfile(output_dir, 'pressure_field.nii'), header);

    focal_vol_bin = int16(p_out > 0.71*max_pressure);
    cc = bwconncomp(focal_vol_bin);
    focal_vol_lcc = int16(zeros(size(p_out)));
    focal_vol_lcc(cc.PixelIdxList{1}) = 1;
    header.Datatype='int16'; header.BitsPerPixel=16;
    niftiwrite(focal_vol_lcc, fullfile(output_dir, 'focal_volume_bin.nii'), header);

    [max_pressure, ~] = max(p_out(logical(focal_vol_lcc))); % [Pa]
    idx = find(p_out==max_pressure);
    [mx, my, mz] = ind2sub(size(p_out), idx);

    disp(['PPW = ' num2str(ppw)])
    disp(['CFL = ' num2str(cfl)])
    disp(['Coordinates of max pressure: [' num2str(mx) ', ' num2str(my) ', ' num2str(mz) ']'])
    disp(['Max Pressure = ' num2str(max_pressure * 1e-6) ' MPa'])
    disp(['MI = ' num2str(MI)])
    disp(['Isppa = ' num2str(Isppa) ' W/cm2'])
    disp(['Pressure at focus = ' num2str(p_focus * 1e-6) ' MPa']);
    disp(['Isppa at focus = ' num2str(isppa_focus) ' W/cm2'])
    disp(['-6dB focal volume = ' num2str(focal_vol) ' mm3'])
    disp(' ')

    clear source sensor_data;
    save(fullfile(output_dir, 'acoustic_sim_output.mat'));

    sim_output_filename = fullfile(output_dir, 'simulation_results.csv');
    if ~exist(sim_output_filename, 'file')
        disp('Result file does not exist, creating file.')
        fileID = fopen(sim_output_filename, 'w' );
        fprintf(fileID, '%s\n', ...
            ['Output directory,Focus coordinates,Bowl coordinates,Focus depth,' ...
            'PPW,CFL,PPP,Coordinates of max pressure,' ...
            'Max Pressure (MPa),MI,Isppa (W/cm2),' ...
            'Pressure at focus (MPa),Isppa at focus (W/cm2),' ...
            '-6dB focal volume (mm3)']);
        fclose(fileID);
    end

    fileID = fopen(sim_output_filename,'a');
    fprintf(fileID,['%s,%d %d %d,%d %d %d,%d,' ...
        '%f,%f,%f,%d %d %d,' ...
        '%f,%f,%f,%f,%f,%f\n'], ...
        output_dir, focus_coords_in, bowl_coords_in, focus_depth,  ...
        ppw, cfl, ppp, mx, my, mz, ...
        max_pressure * 1e-6, MI, Isppa, p_focus * 1e-6, isppa_focus, focal_vol);
    fclose(fileID);
end

end 
