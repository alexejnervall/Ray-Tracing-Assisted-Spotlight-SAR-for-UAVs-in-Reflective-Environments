clearvars; close all; clc

%% %% ----------- RADAR PARAMETERS ----------- %% %%

c = physconst('LightSpeed');

fc = 10e9;                                      % Carrier frequency [Hz]

rangeResolution = 1;                            % [m]
crossRangeResolution = 1;                       % [m]

bw = c/(2*rangeResolution);                     % bandwidth [Hz]

prf = 1500;                                      % [Hz]
aperture = 1;                                   % [m^2]

tpd = 3*10^-6;                                  % Pulse duration [s]
fs = 300*10^6;                                  % Sampling frequency [Hz]
lambda = c/fc;


waveform = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', tpd, 'PRF', prf,...
    'SweepBandwidth', bw);

%% %% ----------- RADAR PLATFORM MOTION ----------- %% %%

speed = 50;
% [m/s]  
maxRange =  700; % even if it should work in theory, the image still comes out very distorted 
synL0 = (lambda*maxRange)/(2*crossRangeResolution)
flightDuration0 = ceil(synL0/speed) % [s]


flightDuration = 1;
synL = flightDuration*speed
prfMin = (2*speed*synL)/(lambda*maxRange)
  

x0 = -0;
y0 = -(synL/2);
z0 = 200;

radarPlatform  = phased.Platform('InitialPosition', [x0; y0; z0], 'Velocity', [0; speed; 0]);

sigma_x = 0.01;   % meters (cross-range disturbance)
sigma_z = 0;   %m
mode = true;     %enable disturbances or not

% Azimuth sampling 
slowTime = 1/prf;
n = flightDuration/slowTime +1                 % Transmitted pulses
slowTimeVec = linspace(0,flightDuration ,n)';   % Slow time vector 

% Range sampling

truncRangesamples = ceil((2*maxRange/c)*fs)
fastTime = (0:1/fs:(truncRangesamples-1)/fs);


%% %% ----------- TARGET DEFINITION ----------- %% %%
enableMapTarget = false;
if enableMapTarget
    scene = uavScenario(ReferenceLocation = [59.402979 17.956008 0]);
    
    xLim = [-100 100];
    yLim = [-100 100];
    
    addMesh(scene,"buildings",{"map2.osm",xLim,yLim,"auto"},[0.6 0.6 0.6]); 
    
    %----------------------------------------------
    points = [];
    
    spacing = 5;   % controls density on surfaces
    
    radarPos = [x0, y0, z0];
    
    for k = 1:numel(scene.Meshes)
    
        verts = scene.Meshes{k}.Vertices;
        faces = scene.Meshes{k}.Faces;
    
        for f = 1:size(faces,1)
    
            v1 = verts(faces(f,1),:);
            v2 = verts(faces(f,2),:);
            v3 = verts(faces(f,3),:);
    
            % --- Compute face normal ---
            n = cross(v2 - v1, v3 - v1);
            if norm(n) == 0
                continue;
            end
            n = n / norm(n);
    
            % --- Face center ---
            center = (v1 + v2 + v3)/3;
    
            % --- Vector toward radar ---
            r = radarPos - center;
    
            % --- Keep ALL faces facing radar (any angle) ---
            if dot(n, r) <= 0
                continue; % backface → skip
            end
    
            % --- Sample triangle surface ---
            % Estimate number of samples based on triangle size
            area = 0.5 * norm(cross(v2 - v1, v3 - v1));
            npts = max(3, ceil(area / (spacing^2)));
    
            for i = 1:npts
                % Random barycentric sampling (uniform over triangle)
                a = rand;
                b = rand;
    
                if a + b > 1
                    a = 1 - a;
                    b = 1 - b;
                end
    
                p = v1 + a*(v2 - v1) + b*(v3 - v1);
                points = [points; p];
            end
    
        end
    end
    
    % Remove duplicates (optional but good)
    points = unique(points,'rows');
    
    targetpos = points';
    targetvel = zeros(size(targetpos));

    % Precompute radar trajectory
    radarPosHist = zeros(3, n);
    
    radarPlatformTemp = phased.Platform( ...
        'InitialPosition', [x0; y0; z0], ...
        'Velocity', [0; speed; 0]);
    
    for ii = 1:n
        [pos, ~] = radarPlatformTemp(slowTime);
        radarPosHist(:,ii) = pos;
    end

    % ----------- GROUND TRUTH VISUALIZATION ----------- %
    figure
    h = axes;
    hold on
    grid on
    axis equal
    
    % Targets in 3D 
    scatter3(targetPos(2,:), targetPos(1,:), targetPos(3,:), 50, 'filled', 'MarkerFaceColor','r');
    
    plot3(radarPosHist(2,:), radarPosHist(1,:), radarPosHist(3,:), 'b', 'LineWidth', 2)
    
    margin = 10;  
    x_min = min(targetPos(2,:)) - margin;
    x_max = max(targetPos(2,:)) + margin;
    y_min = min(targetPos(1,:)) - margin;
    y_max = max(targetPos(1,:)) + margin;
    z_min = min(targetPos(3,:)) - margin;
    z_max = max(targetPos(3,:)) + margin;
    
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    zlim([z_min, z_max])
    
    xlabel('Cross-Range (y)')
    ylabel('Range (x)')
    zlabel('Height (z)')
    title('3D Ground Truth')
    set(h,'YDir','reverse')
    
    view(45,30)  
end

%% ----------------------------------------------

% TESTING WITH A CUBE FOR SPEED 

% ----------------------------------------------
enableCubeTarget = false;
if enableCubeTarget
    points = [];
    
    spacing = 3;   % controls point density along edges
    
    cubeSize = 16;    
    halfSize = cubeSize / 2;
    
    % --- Define cube vertices ---
    verts = [
        -halfSize, -halfSize, -halfSize;
         halfSize, -halfSize, -halfSize;
         halfSize,  halfSize, -halfSize;
        -halfSize,  halfSize, -halfSize;
        -halfSize, -halfSize,  halfSize;
         halfSize, -halfSize,  halfSize;
         halfSize,  halfSize,  halfSize;
        -halfSize,  halfSize,  halfSize
    ];
    
    theta = deg2rad(60);  % 45 degrees
    
    Rz = [
        cos(theta), -sin(theta), 0;
        sin(theta),  cos(theta), 0;
        0,           0,          1
    ];
    
    % Rotate all vertices
    verts = (Rz * verts')';
    
    % --- Define cube edges (pairs of vertex indices) ---
    edges = [
        1 2; 2 3; 3 4; 4 1;  % bottom square
        5 6; 6 7; 7 8; 8 5;  % top square
        1 5; 2 6; 3 7; 4 8   % vertical edges
    ];
    
    for e = 1:size(edges,1)
        v1 = verts(edges(e,1),:);
        v2 = verts(edges(e,2),:);
    
        % --- Compute number of points along edge ---
        edgeLength = norm(v2 - v1);
        npts = max(2, ceil(edgeLength / spacing));
    
        for i = 0:npts
            t = i / npts;
            p = v1 + t*(v2 - v1);
            points = [points; p];
        end
    end
    
    % Remove duplicates
    points = unique(points,'rows');
    
    targetPos = points';
    targetVel = zeros(size(targetPos));

    % ----------- GROUND TRUTH VISUALIZATION ----------- %
    % Precompute radar trajectory
    radarPosHist = zeros(3, n);
    
    radarPlatformTemp = phased.Platform( ...
        'InitialPosition', [x0; y0; z0], ...
        'Velocity', [0; speed; 0]);
    
    for ii = 1:n
        [pos, ~] = radarPlatformTemp(slowTime);
        radarPosHist(:,ii) = pos;
    end

    
    figure
    h = axes;
    hold on
    grid on
    axis equal
    
    % Targets in 3D 
    scatter3(targetPos(2,:), targetPos(1,:), targetPos(3,:), 50, 'filled', 'MarkerFaceColor','r');
    
    plot3(radarPosHist(2,:), radarPosHist(1,:), radarPosHist(3,:), 'b', 'LineWidth', 2)
    
    margin = 10;  
    x_min = min(targetPos(2,:)) - margin;
    x_max = max(targetPos(2,:)) + margin;
    y_min = min(targetPos(1,:)) - margin;
    y_max = max(targetPos(1,:)) + margin;
    z_min = min(targetPos(3,:)) - margin;
    z_max = max(targetPos(3,:)) + margin;
    
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    zlim([z_min, z_max])
    
    xlabel('Cross-Range (y)')
    ylabel('Range (x)')
    zlabel('Height (z)')
    title('3D Ground Truth')
    set(h,'YDir','reverse')
    
    view(45,30)  
end
% ----------------------------------------------

% END OF CUBE TESTING 

% ----------------------------------------------


%% TESTING WITH 3 point target
enable3Target = true;
if enable3Target
    targetPos= [100,0,0;110,-30,0; 130, -15, 0]';
    targetVel = zeros(size(targetPos));

    % ----------- GROUND TRUTH VISUALIZATION ----------- %
    figure
    radarPosHist = zeros(3, n);
    
    radarPlatformTemp = phased.Platform( ...
        'InitialPosition', [x0; y0; z0], ...
        'Velocity', [0; speed; 0]);
    
    for ii = 1:n
        [pos, ~] = radarPlatformTemp(slowTime);
        radarPosHist(:,ii) = pos;
    end

    
    figure
    h = axes;
    hold on
    grid on
    axis equal
    
    % Targets in 3D 
    scatter3(targetPos(2,:), targetPos(1,:), targetPos(3,:), 50, 'filled', 'MarkerFaceColor','r');
    
    plot3(radarPosHist(2,:), radarPosHist(1,:), radarPosHist(3,:), 'b', 'LineWidth', 2)
    
    margin = 10;  
    x_min = min(targetPos(2,:)) - margin;
    x_max = max(targetPos(2,:)) + margin;
    y_min = min(targetPos(1,:)) - margin;
    y_max = max(targetPos(1,:)) + margin;
    z_min = min(targetPos(3,:)) - margin;
    z_max = max(targetPos(3,:)) + margin;
    
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    zlim([z_min, z_max])
    
    xlabel('Cross-Range (y)')
    ylabel('Range (x)')
    zlabel('Height (z)')
    title('3D Ground Truth')
    set(h,'YDir','reverse')
    
    view(45,30)

end
%% --------- Geometry Calculations ---------%%

sceneCenter = mean(targetPos, 2);
lookVec = sceneCenter - radarPosHist;
lookVecNorm = lookVec ./ vecnorm(lookVec);
Rc = norm(sceneCenter - radarPosHist(:,round(end/2)))
xCenter = mean(targetPos(1,:))
yCenter = mean(targetPos(2,:))
AoA = atand(z0/xCenter)
squintAngle = 0;
%% %% ----------- SAR SIGNAL SIMULATION USING RAY TRACING ----------- %% %%

% Target object applies amplitude scaling and phase shift to reflected signal 
pointTargets = phased.Platform('InitialPosition', targetPos,'Velocity',targetVel);

pm = propagationModel("raytracing", ...
    Method="sbr", ...
    CoordinateSystem="cartesian", ...
    MaxNumReflections=1, ...
    MaxNumDiffractions=0);

rxSig = zeros(truncRangesamples,n);

%Initilazing for disturbances 

radarPosNominal = zeros(3,n);

for ii = 1:n
    [pos, ~] = radarPlatform(slowTime);
    radarPosNominal(:,ii) = pos;
end

[radarPosDisturbed, trajError] = generateDisturbedTrajectory( ...
    radarPosNominal, sigma_x, sigma_z, mode);

for ii = 1:n
    
    fprintf('Pulse %d / %d\n', ii, n)

    sig = waveform();
    sig = sig(1:truncRangesamples);


    radarpos = radarPosDisturbed(:,ii);
    [targetPos,targetVel] = pointTargets(slowTime);
 

    tx = txsite("cartesian", "AntennaPosition", radarpos(:), "TransmitterFrequency", fc);

    for k = 1:size(targetPos,2)

        rx = rxsite("cartesian", ...
            "AntennaPosition", targetPos(:,k));

        rays = raytrace(tx,rx,pm);

        if isempty(rays)
            continue
        end

        rayset = rays{1};
        
        for r = 1:length(rayset)

            dist = rayset(r).PropagationDistance;
            tau = 2*dist/c;
            delay = round(tau*fs);
            loss = db2pow(-rayset(r).PathLoss);
            phase = loss * exp(-1i*4*pi*fc*dist/c);
            
            % Add delayed echo to received signal
            if delay < truncRangesamples
            
                valid = 1:(truncRangesamples-delay);
                rxSig(delay + valid, ii) = rxSig(delay + valid, ii) + phase*sig(valid);
            
            end
        end
    end
end


%% %% ----------- SIGNAL PROCESSING ----------- %% %%

% Visualization of raw SAR data

figure(2)
imagesc(real(rxSig));title('SAR Raw Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')
%colormap('gray')

% Range compression

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
matchingCoeff = getMatchedFilter(waveform);
[cdata, rnggrid] = pulseCompression(rxSig, matchingCoeff);

% Range compressed data
figure (3);
imagesc(real(cdata));
title('SAR Range Compressed Data ');
xlabel('Cross-Range Samples');
ylabel('Range Samples');
%colormap('gray')

%% --------------- Error Plots -------------------%%
figure
plot(trajError(1,:))
title('Cross-range disturbance')

figure
plot(trajError(3,:))
title('Altitude disturbance')
figure(2);
%% %% ----------- IMAGE FORMATION ----------- %% %%

% Computing reference range 
Rc = norm(sceneCenter - radarPosHist(:,round(end/2)));

% Azimuth compression
rma_processed = RMA(cdata,fastTime,fc,fs,prf,speed,n,c,Rc,squintAngle);

% Final image, compressed in range & azimuth
figure(4)
imagesc(abs(rma_processed).')
title('Full SAR Image')
xlabel('Cross-Range Samples')
ylabel('Range Samples')
%colormap('gray')