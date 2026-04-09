clearvars; close all; clc

%% %% ----------- RADAR PARAMETERS ----------- %% %%

c = physconst('LightSpeed');

fc = 10e9;                          % Carrier frequency [Hz]

dR = 1;                             % [m]
dA = 1;                             % [m]

bw = c/(2*dR);                      % bandwidth [Hz]

prf = 1500;                         % [Hz]
aperture = 0.673;                   % [m^2]

tpd = 3*10^-6;                      % Pulse duration [s]
fs = 300*10^6;                      % Sampling frequency [Hz]
lambda = c/fc;                      % Wavelength [m]


waveform = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', tpd, 'PRF', prf,...
    'SweepBandwidth', bw);

%% %% ----------- RADAR PLATFORM MOTION ----------- %% %%

v = 50;                                     
Rmax = 700; 
L0syn = (lambda*Rmax)/(2*dA);
T0 = ceil(L0syn/v);                           


T = 1;
Lsyn = T*v;
prfMin = (2*v*Lsyn)/(lambda*Rmax);

x0 = 0;
y0 = -(Lsyn/2);
z0 = 200;

radarPlatform  = phased.Platform('InitialPosition', [x0; y0; z0], 'Velocity', [0; v; 0]);

% Azimuth sampling 
slowTime = 1/prf;
n = T/slowTime +1;                  % Transmitted pulses
slowTimeVec = linspace(0,T , n)';   % Slow time vector 

% Range sampling
truncRangesamples = ceil((2*Rmax/c)*fs);
fastTime = (0:1/fs:(truncRangesamples-1)/fs);


%% %% ----------- TARGET DEFINITION ----------- %% %%

mapTarget = false;
cubeTarget = true;

points = pointCloudGeneration(x0, y0, z0, mapTarget, cubeTarget);

targetPos = points';
targetVel = zeros(size(targetPos));

% Target object applies amplitude scaling and phase shift to reflected signal 
target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', ones(1,size(targetPos,2)));
pointTargets = phased.Platform('InitialPosition', targetPos,'Velocity',targetVel);


%% %% ----------- GROUND TRUTH VISUALIZATION ----------- %% %%

scatter3(targetPos(1,:), targetPos(2,:), targetPos(3,:), 20, 'filled');
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Point Cloud of Cube Edges');


% Precompute radar trajectory
radarPosHist = zeros(3, n);

radarPlatformTemp = phased.Platform( ...
    'InitialPosition', [x0; y0; z0], ...
    'Velocity', [0; v; 0]);

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

%% %% ----------- GEOMETRY CALCULATIONS ----------- %% %%


sceneCenter = mean(targetPos, 2);
lookVec = sceneCenter - radarPosHist;
lookVecNorm = lookVec ./ vecnorm(lookVec);
Rc = norm(sceneCenter - radarPosHist(:, round(end/2)));
xc = mean(targetPos(1, :));
yc = mean(targetPos(2, :));
AoA = atand(z0/xc);
squintAngle = 0;

%% %% ----------- SAR SIGNAL SIMULATION USING RAY TRACING ----------- %% %%

% Disturbances 

sx = 0.01;
sz = 0.0; 
disturbances = false;

radarPosNom = zeros(3, n);

for jj = 1:n
    [pos, ~] = radarPlatform(slowTime);
    radarPosNom(:, jj) = pos;
end

[radarPosDisturbed, trajError] = trajectoryGeneration(radarPosNom, sx, sz, disturbances);

% --- 

pm = propagationModel("raytracing", ...
    Method="sbr", ...
    CoordinateSystem="cartesian", ...
    MaxNumReflections=1, ...
    MaxNumDiffractions=0);

rxSig = zeros(truncRangesamples,n);

for ii = 1:n
    
    fprintf('Pulse %d / %d\n', ii, n)

    sig = waveform();
    sig = sig(1:truncRangesamples);

    radarPos = radarPosDisturbed(:, ii);
    %[radarpos, radarvel] = radarPlatform(slowTime);
    [targetPos,targetVel] = pointTargets(slowTime);

    tx = txsite("cartesian", "AntennaPosition", radarPos(:), "TransmitterFrequency", fc);

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
colormap('gray')

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
colormap('gray')


%% %% ----------- IMAGE FORMATION ----------- %% %%

% Azimuth compression
rmaProcessed = RMA(cdata, fastTime, fc, fs, prf, v, n, c, Rc, squintAngle);

sarImage = abs(rmaProcessed);

% Final image, compressed in range & azimuth
figure(4)
imagesc(sarImage.')
title('Full SAR Image')
xlabel('Cross-Range Samples')
ylabel('Range Samples')
colormap('gray')

%% %% ----------- IMAGE ANALYSIS ----------- %% %%

rmsContrast = computeRMSContrast(sarImage)

entropy = entropy(sarImage)


%% %% ----------- ERROR PLOTS ----------- %% %%

figure(5)
plot(trajError(1, :))
title('Crossrange Disturbance')
xlabel('Pulse Number')         
ylabel('Error [m]')             
legend('Crossrange Error')      

figure(6)
plot(trajError(3, :))
title('Altitude Disturbance')
xlabel('Pulse Number')         
ylabel('Error [m]')             
legend('Altitude Error')     