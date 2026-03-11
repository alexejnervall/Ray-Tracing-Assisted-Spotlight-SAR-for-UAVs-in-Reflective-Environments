clearvars; close all; clc

%% Radar params

c = physconst('LightSpeed');

fc = 4e9;                                       % Carrier frequency [Hz]

rangeResolution = 3;                            % [meters]
crossRangeResolution = 3;                       % [meters]

bw = c/(2*rangeResolution);                     % bandwidth

prf = 1000;                                     % [Hz]
aperture = 4;                                   % [sq. meters]
tpd = 3*10^-6;                                  % Pulse duration [s]
fs = 120*10^6;                                  % Sampling frequency [Hz]


waveform = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', tpd, 'PRF', prf,...
    'SweepBandwidth', bw);

%% Radar platform motion 

speed = 100;                                    % [m/s]  
flightDuration = 4;                             % [s]
radarPlatform  = phased.Platform('InitialPosition', [0;-600;500], 'Velocity', [0; speed; 0]);
slowTime = 1/prf;
numpulses = flightDuration/slowTime +1;         % Transmitted pulses
eta1 = linspace(0,flightDuration ,numpulses)';  % Slow time vector 

% Range sampling
maxRange = 2500;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime = (0:1/fs:(truncrangesamples-1)/fs);

% Set the reference range for the cross-range processing.
Rc = 1e3;                                       % [m]

%% Tx and Rx config 

% Creating single antenna element with cosine radiation pattern 
antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);

% Converting physical antenna area to antenna gain in dB 
antennaGain = aperture2gain(aperture,c/fc); 

% Transmitter: scaling signal amplitude by power and antenna gain
transmitter = phased.Transmitter('PeakPower', 1e3, 'Gain', antennaGain);

% Radiator: converting  digital waveform to radiated EM field in 3D space
radiator = phased.Radiator('Sensor', antenna,'OperatingFrequency', fc, 'PropagationSpeed', c);

% Collector: converting incoming EM waves from targets into a baseband signal at the receiver
collector = phased.Collector('Sensor', antenna, 'PropagationSpeed', c,'OperatingFrequency', fc);

% Receiver: taking collected signal and outputs digitized radar echoes 
receiver = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);

%% Target

% Stationary targets 
targetpos= [900,0,0;1000,-30,0; 850, -15, 0]';
targetvel = [0,0,0;0,0,0; 0,0,0]';

squintangle = atand(600/950);

% Target object applies amplitude scaling and phase shift to reflected signal 
target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', ones(1,size(targetpos,2)));
pointTargets = phased.Platform('InitialPosition', targetpos,'Velocity',targetvel);

% Ground truth plot (plotting all targets regardless of number) 

figure(1)
h = axes;

plot(targetpos(2,:), targetpos(1,:), '*')

set(h,'Ydir','reverse')

xlim([-50 10])
ylim([800 1200])   

title('Ground Truth')
ylabel('Range')
xlabel('Cross-Range')

% Signal simulation using ray tracing

pm = propagationModel("raytracing", ...
    Method="sbr", ...
    CoordinateSystem="cartesian", ...
    MaxNumReflections=1, ...
    MaxNumDiffractions=0);

% Init 2D matrix for SAR data 
rxsig = zeros(truncrangesamples,numpulses);

% Looping over all radar pulses 
for ii = 1:numpulses

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

            % Generating an LFM pulse 
            sig = waveform();
            sig = sig(1:truncrangesamples);
            
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


kc = (2*pi*fc)/c;
 
% Phase correction for Doppler due to squint angle 
rxsig=rxsig.*exp(-1i.*2*(kc)*sin(deg2rad(squintangle))*repmat(speed*eta1,1,truncrangesamples)).';

% Visualization 

figure(2)
imagesc(real(rxsig));title('SAR Raw Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

% Range compression

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
matchingCoeff = getMatchedFilter(waveform);
[cdata, rnggrid] = pulseCompression(rxsig, matchingCoeff);


figure (3);
imagesc(real(cdata));
title('SAR Range Compressed Data ');
xlabel('Cross-Range Samples');
ylabel('Range Samples');

% Azimuth compression
rma_processed = helperSquintRangeMigration(cdata,fastTime,fc,fs,prf,speed,numpulses,c,Rc,squintangle);

% Final image (Plotting whole area so we don't exclude any targets)

figure(4)
imagesc(abs(rma_processed).')
title('Full SAR Image')
xlabel('Cross-Range Samples')
ylabel('Range Samples')


%% Azimuth compression helper function 

function azcompresseddata = helperSquintRangeMigration(sigData,fastTime,fc,fs,prf,speed,numPulses,c,Rc,squintangle)

% Creating RF frequency grid for sampled signal 
frequencyRange = linspace(fc-fs/2,fc+fs/2,length(fastTime));

% Range spacial frequency axis 
krange = 2*(2*pi*frequencyRange)/c;

% Azimuth spacial frequency axis 
kaz = 2*pi*linspace(-prf/2,prf/2,numPulses)./speed;

% Carrier spacial frequency 
kc = 2*pi*fc/3e8;
kazimuth = kaz.';

% Doppler center shift 
kus=2*(kc)*sin(deg2rad(squintangle));

% SAR geometry equation
kx = krange.^2-(kazimuth+kus).^2;

% Converting squint angle to rad
thetaRc = deg2rad(squintangle);
kx = sqrt(kx.*(kx > 0));

% Phase correction 
kFinal = exp(1i*(kx.*cos(thetaRc)+(kazimuth).*sin(thetaRc)).*Rc);
kfin = kx.*cos(thetaRc)+(kazimuth+kus).*sin(thetaRc);

% Convert to 2D frequency domain (2D FFT) 
sdata =fftshift(fft(fftshift(fft(sigData,[],1),1),[],2),2);

% Applying phase compensation 
fsmPol = (sdata.').*kFinal;

% Stolt interpolation 
stoltPol = fsmPol;
for i = 1:size((fsmPol),1)
    stoltPol(i,:) = interp1(kfin(i,:),fsmPol(i,:),krange(1,:));
end
stoltPol(isnan(stoltPol)) = 1e-30;

% 2D IFFT 
azcompresseddata = ifftshift(ifft2(stoltPol),2);

end