clearvars; close all; clc

%% Radar params

c = physconst('LightSpeed');

fc = 10e9;                                       % Carrier frequency [Hz]

rangeResolution = 1;                            % [meters]
crossRangeResolution = 1;                       % [meters]

bw = c/(2*rangeResolution)                     % bandwidth

prf = 500;                                     % [Hz]
aperture = 1;                                   % [sq. meters]

tpd = 3*10^-6;                                  % Pulse duration [s]
fs = 300*10^6;                                  % Sampling frequency [Hz]


waveform = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', tpd, 'PRF', prf,...
    'SweepBandwidth', bw);

%% Radar platform motion 

speed = 30;                                    % [m/s]  
flightDuration = 4;                             % [s]
radarPlatform  = phased.Platform('InitialPosition', [400;-80;30], 'Velocity', [0; speed; 0]);
slowTime = 1/prf;
numpulses = flightDuration/slowTime +1;         % Transmitted pulses
eta1 = linspace(0,flightDuration ,numpulses)';  % Slow time vector 

% Range sampling
maxRange = 2000;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime = (0:1/fs:(truncrangesamples-1)/fs);

%% Target

% Stationary targets 
targetpos= [900,0,0;800,-30,0; 850, -15, 0]';
targetvel = [0,0,0;0,0,0; 0,0,0]';

%squintangle = atand(235/950);
squintangle = 0;

% Target object applies amplitude scaling and phase shift to reflected signal 
target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', ones(1,size(targetpos,2)));
pointTargets = phased.Platform('InitialPosition', targetpos,'Velocity',targetvel);

%% Ground truth plot 

% Precompute radar trajectory
radarPosHist = zeros(3, numpulses);

radarPlatformTemp = phased.Platform( ...
    'InitialPosition', [400;-80;30], ...
    'Velocity', [0; speed; 0]);

for ii = 1:numpulses
    [pos, ~] = radarPlatformTemp(slowTime);
    radarPosHist(:,ii) = pos;
end

sceneCenter = mean(targetpos, 2)
lookVec = sceneCenter - radarPosHist;
lookVecNorm = lookVec ./ vecnorm(lookVec);
%squintangle = atan2d(lookVec(2), lookVec(1))


figure(1)
h = axes;

plot(targetpos(2,:), targetpos(1,:), '*', 'MarkerSize', 8)
set(h,'YDir','reverse')

% Automatically compute axis limits with margin
margin = 10;  % meters
x_min = min(targetpos(2,:)) - margin;
x_max = max(targetpos(2,:)) + margin;
y_min = min(targetpos(1,:)) - margin;
y_max = max(targetpos(1,:)) + margin;

xlim([x_min, x_max])
ylim([y_min, y_max])

title('Ground Truth')
ylabel('Range')
xlabel('Cross-Range')
grid on
axis equal
%%
figure(1)
hold on
grid on
axis equal

% --- Targets ---
scatter3(targetpos(2,:), targetpos(1,:), targetpos(3,:), 50, 'filled')

% --- Radar trajectory ---
plot3(radarPosHist(2,:), radarPosHist(1,:), radarPosHist(3,:), ...
    'b', 'LineWidth', 2)

% --- Radar positions (sparse markers) ---
idx = 1:round(numpulses/20):numpulses; % reduce clutter
scatter3(radarPosHist(2,idx), radarPosHist(1,idx), radarPosHist(3,idx), ...
    20, 'b', 'filled')

% --- Look direction vectors (quiver) ---
scale = 200; % adjust for visibility

quiver3( ...
    radarPosHist(2,idx), ...
    radarPosHist(1,idx), ...
    radarPosHist(3,idx), ...
    lookVecNorm(2,idx)*scale, ...
    lookVecNorm(1,idx)*scale, ...
    lookVecNorm(3,idx)*scale, ...
    0, 'r')

xlabel('Cross-range (y)')
ylabel('Range (x)')
zlabel('Height (z)')

title('Ground Truth with Radar Trajectory and Look Direction')

legend('Targets', 'Radar Path', 'Radar Samples', 'Look Direction')
view(45,30)


% Set the reference range for the cross-range processing.
%Rc = 1e3;                                       % [m]

sceneCenter = mean(targetpos,2);
Rc = norm(sceneCenter - radarPosHist(:,round(end/2)))
%% Tx and Rx config 

% Creating single antenna element with cosine radiation pattern 
antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);

% Converting physical antenna area to antenna gain in dB 
antennaGain = aperture2gain(aperture,c/fc); 

% Transmitter: scaling signal amplitude by power and antenna gain
transmitter = phased.Transmitter('PeakPower', 1e3, 'Gain', antennaGain);


%% Signal simulation using ray tracing

pm = propagationModel("raytracing", ...
    Method="sbr", ...
    CoordinateSystem="cartesian", ...
    MaxNumReflections=1, ...
    MaxNumDiffractions=0);

% Init 2D matrix for SAR data 
rxsig = zeros(truncrangesamples,numpulses);

% Looping over all radar pulses 
for ii = 1:numpulses
    
    fprintf('Pulse %d / %d\n', ii, numpulses)

    % Generating an LFM pulse 
    sig = waveform();
    sig = sig(1:truncrangesamples);
    sig = transmitter(sig);

    % Update radar platform and target position
    [radarpos, radarvel] = radarPlatform(slowTime);
    [targetpos,targetvel] = pointTargets(slowTime);

    % Modeling the radar as Tx in 3D space 
    tx = txsite("cartesian", ...
        "AntennaPosition", radarpos(:), ...
        "TransmitterFrequency", fc);

    % Looping over each target
    for k = 1:size(targetpos,2)

        % Target as receiver of the pulse 
        rx = rxsite("cartesian", ...
            "AntennaPosition", targetpos(:,k));

        % Ray tracing
        rays = raytrace(tx,rx,pm);

        if isempty(rays)
            continue
        end

        rayset = rays{1};
        
        % Loop over rays for each target
        for r = 1:length(rayset)

            % Propagation distance along the ray
            dist = rayset(r).PropagationDistance;

            % Round trip delay 
            tau = 2*dist/c;
            
            % Computing delay and phase 
            delay = round(tau*fs);
            loss = db2pow(-rayset(r).PathLoss);
            phase = loss * exp(-1i*4*pi*fc*dist/c);
            
            % Add delayed echo to raw data 
            if delay < truncrangesamples
            
                valid = 1:(truncrangesamples-delay);
            
                rxsig(delay + valid, ii) = ...
                    rxsig(delay + valid, ii) + phase*sig(valid);
            
            end
        end
    end
end

%%
kc = (2*pi*fc)/c;
 
% Phase correction for Doppler due to squint angle 
rxsig=rxsig.*exp(-1i.*2*(kc)*sin(deg2rad(squintangle))*repmat(speed*eta1,1,truncrangesamples)).';

% Visualization of raw SAR data

figure(2)
imagesc(real(rxsig));title('SAR Raw Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')
colormap('gray')

% Range compression

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
matchingCoeff = getMatchedFilter(waveform);
[cdata, rnggrid] = pulseCompression(rxsig, matchingCoeff);


figure (3);
imagesc(real(cdata));
title('SAR Range Compressed Data ');
xlabel('Cross-Range Samples');
ylabel('Range Samples');
colormap('gray')

% Azimuth compression
rma_processed = RMA(cdata,fastTime,fc,fs,prf,speed,numpulses,c,Rc,squintangle);

% Final image

figure(4)
imagesc(abs(rma_processed).')
title('Full SAR Image')
xlabel('Cross-Range Samples')
ylabel('Range Samples')
colormap('gray')