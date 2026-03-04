clearvars; close all; clc

%% ============================================================
% Load OSM Map (Geographic)
%% ============================================================
viewer = siteviewer(Buildings="map.osm");

pm = propagationModel("raytracing", ...
    Method="sbr", ...
    MaxNumReflections=1, ...
    MaxNumDiffractions=0);

%% ============================================================
% Radar Parameters
%% ============================================================
c  = physconst('LightSpeed');
fc = 10e9;

rangeResolution = 3;
bw = c/(2*rangeResolution);

prf = 400;               % MUST divide fs
tpd = 3e-6;
fs  = 40e6;              % 40e6 / 400 = 100000 ✔ integer

waveform = phased.LinearFMWaveform( ...
    'SampleRate',fs, ...
    'PulseWidth',tpd, ...
    'PRF',prf, ...
    'SweepBandwidth',bw);

speed = 100;             
flightDuration = 2;      % keep short for speed
numpulses = floor(flightDuration * prf);
eta = (0:numpulses-1)'/prf;

maxRange = 2500;
truncrangesamples = ceil((2*maxRange/c)*fs);

%% ============================================================
% Reference Geographic Location
%% ============================================================
lat0 = 59.405023;
lon0 = 17.947607;
platformHeight = 500;     
earthRadius = 6378137;

%% ============================================================
% SAR Data Collection (FAST RAY-DRIVEN)
%% ============================================================
rxsig = zeros(truncrangesamples,numpulses);

disp("Starting SAR collection (fast version)...")

for ii = 1:numpulses

    % Platform motion northbound
    dNorth = speed * eta(ii);
    radarLat = lat0 + (dNorth/earthRadius)*(180/pi);
    radarLon = lon0;

    tx = txsite( ...
        Latitude=radarLat, ...
        Longitude=radarLon, ...
        AntennaHeight=platformHeight, ...
        TransmitterFrequency=fc);

    % Raytrace only to scene center
    tgt = rxsite( ...
        Latitude=lat0, ...
        Longitude=lon0, ...
        AntennaHeight=0);

    rays = raytrace(tx,tgt,pm);

    sig = waveform();
    sig = sig(1:truncrangesamples);

    rxPulse = zeros(truncrangesamples,1);

    if ~isempty(rays) && ~isempty(rays{1})

        for r = 1:length(rays{1})

            tau = rays{1}(r).PropagationDelay;
            L   = rays{1}(r).PathLoss;
            phi = rays{1}(r).PhaseShift;

            % Two-way propagation
            tau = 2*tau;
            L   = 2*L;

            gain = 10^(-L/20);
            sampleDelay = round(tau * fs);

            if sampleDelay < truncrangesamples
                rxPulse(sampleDelay+1:end) = ...
                    rxPulse(sampleDelay+1:end) + ...
                    gain * sig(1:end-sampleDelay) ...
                    * exp(1j*phi);
            end
        end
    end

    rxsig(:,ii) = rxPulse;

    if mod(ii,50)==0
        fprintf("Pulse %d / %d\n",ii,numpulses);
    end
end

disp("SAR collection complete.")

%% ============================================================
% Squint Compensation
%% ============================================================
squintangle = 33;
kc = (2*pi*fc)/c;

dopplerPhase = exp( ...
    -1i * 2 * kc * sin(deg2rad(squintangle)) ...
    * speed * eta.' );     % 1 × numPulses

rxsig = rxsig .* dopplerPhase;   % auto-expands over rows
%% ============================================================
% Range Compression
%% ============================================================
pulseCompression = phased.RangeResponse( ...
    'RangeMethod','Matched filter', ...
    'PropagationSpeed',c, ...
    'SampleRate',fs);

matchingCoeff = getMatchedFilter(waveform);
[cdata,~] = pulseCompression(rxsig,matchingCoeff);

figure
imagesc(real(cdata))
title('Range Compressed Data')
xlabel('Cross-Range')
ylabel('Range')
hold on

%% ============================================================
% Range Migration Algorithm
%% ============================================================
Rc = 1000;
fastTime = (0:truncrangesamples-1)/fs;

%rma_processed = helperSquintRangeMigration( ...
    %cdata,fastTime,fc,fs,prf,speed,numpulses,c,Rc,squintangle);

rma_processed = rangeMigrationLFM(cdata, waveform, fc, speed, maxRange);

figure
imagesc(abs(rma_processed.'))
title('Focused Urban SAR Image')
xlabel('Cross-Range')
ylabel('Range')
colormap hot
colorbar

%% ============================================================
% Backprojection Algorithm
%% ============================================================

% The backprojection function doesn't seem to exist anymore. Will have to
% be written manually 

% bpa_processed = helperBackProjection