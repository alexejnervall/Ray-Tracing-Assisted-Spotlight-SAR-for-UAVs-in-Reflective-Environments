clearvars; close all; clc

%% %% ----------- RADAR PARAMETERS ----------- %% %%

c = physconst('LightSpeed');

fc = 10e9;                                      % Carrier frequency [Hz]

rangeResolution = 1;                            % [m]
crossRangeResolution = 1;                       % [m]

bw = c/(2*rangeResolution);                     % bandwidth [Hz]

prf = 500;                                      % [Hz]
aperture = 1;                                   % [m^2]

tpd = 3*10^-6;                                  % Pulse duration [s]
fs = 300*10^6;                                  % Sampling frequency [Hz]


waveform = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', tpd, 'PRF', prf,...
    'SweepBandwidth', bw);

%% %% ----------- RADAR PLATFORM MOTION ----------- %% %%

speed = 100;                                     % [m/s]  
flightDuration = 4;                             % [s]

x0 = 400;
y0 = -200;
z0 = 200;

radarPlatform  = phased.Platform('InitialPosition', [x0; y0; z0], 'Velocity', [0; speed; 0]);

% Azimuth sampling 
slowTime = 1/prf;
n = flightDuration/slowTime +1;                 % Transmitted pulses
slowTimeVec = linspace(0,flightDuration ,n)';   % Slow time vector 

% Range sampling
maxRange = 2000;
truncRangesamples = ceil((2*maxRange/c)*fs);
fastTime = (0:1/fs:(truncRangesamples-1)/fs);

%% %% ----------- TX AND RX CONFIG ----------- %% %%


% Creating single antenna element with cosine radiation pattern 
antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]); % THIS FREQUENCY RANGE IS WRONG NO???!!!!

% Converting physical antenna area to antenna gain in dB 
antennaGain = aperture2gain(aperture,c/fc); 

% Transmitter: scaling signal amplitude by power and antenna gain
transmitter = phased.Transmitter('PeakPower', 1e3, 'Gain', antennaGain);


%% %% ----------- TARGET DEFINITION ----------- %% %%


scene = uavScenario(ReferenceLocation = [59.402979 17.956008 0]);

xLim = [-100 100];
yLim = [-100 100];

addMesh(scene,"buildings",{"map2.osm",xLim,yLim,"auto"},[0.6 0.6 0.6]); 

% ----------------------------------------------
% points = [];
% 
% spacing = 5;   % controls density on surfaces
% 
% radarPos = [x0, y0, z0];
% 
% for k = 1:numel(scene.Meshes)
% 
%     verts = scene.Meshes{k}.Vertices;
%     faces = scene.Meshes{k}.Faces;
% 
%     for f = 1:size(faces,1)
% 
%         v1 = verts(faces(f,1),:);
%         v2 = verts(faces(f,2),:);
%         v3 = verts(faces(f,3),:);
% 
%         % --- Compute face normal ---
%         n = cross(v2 - v1, v3 - v1);
%         if norm(n) == 0
%             continue;
%         end
%         n = n / norm(n);
% 
%         % --- Face center ---
%         center = (v1 + v2 + v3)/3;
% 
%         % --- Vector toward radar ---
%         r = radarPos - center;
% 
%         % --- Keep ALL faces facing radar (any angle) ---
%         if dot(n, r) <= 0
%             continue; % backface → skip
%         end
% 
%         % --- Sample triangle surface ---
%         % Estimate number of samples based on triangle size
%         area = 0.5 * norm(cross(v2 - v1, v3 - v1));
%         npts = max(3, ceil(area / (spacing^2)));
% 
%         for i = 1:npts
%             % Random barycentric sampling (uniform over triangle)
%             a = rand;
%             b = rand;
% 
%             if a + b > 1
%                 a = 1 - a;
%                 b = 1 - b;
%             end
% 
%             p = v1 + a*(v2 - v1) + b*(v3 - v1);
%             points = [points; p];
%         end
% 
%     end
% end
% 
% % Remove duplicates (optional but good)
% points = unique(points,'rows');
% 
% targetpos = points';
% targetvel = zeros(size(targetpos));

% ----------------------------------------------

% TESTING WITH A CUBE FOR SPEED 

% ----------------------------------------------

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

theta = deg2rad(45);  % 45 degrees

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

% ----------------------------------------------

% END OF CUBE TESTING 

% ----------------------------------------------

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
    'Velocity', [0; speed; 0]);

for ii = 1:n
    [pos, ~] = radarPlatformTemp(slowTime);
    radarPosHist(:,ii) = pos;
end

sceneCenter = mean(targetPos, 2);
lookVec = sceneCenter - radarPosHist;
lookVecNorm = lookVec ./ vecnorm(lookVec);
squintAngle = 0;

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


%% %% ----------- SAR SIGNAL SIMULATION USING RAY TRACING ----------- %% %%

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
    sig = transmitter(sig);

    [radarpos, radarvel] = radarPlatform(slowTime);
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

% kc = (2*pi*fc)/c;  % !!!!!!!!!!!!!!!! NOT NEEDED FOR SQUINT = 0, REMOVE AFTER DOUBLE CHECK!!!!!!!!!!!!!!!!

% % Phase correction for Doppler due to squint angle 
% rxSig=rxSig.*exp(-1i.*2*(kc)*sin(deg2rad(squintAngle))*repmat(speed*slowTimeVec,1,truncRangesamples)).';

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
colormap('gray')