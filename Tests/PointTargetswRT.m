clearvars; close all; clc

% Radar params

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

% Radar platform motion 

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

% Tx and Rx config 

antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(aperture,c/fc); 

transmitter = phased.Transmitter('PeakPower', 1e3, 'Gain', antennaGain);
radiator = phased.Radiator('Sensor', antenna,'OperatingFrequency', fc, 'PropagationSpeed', c);

collector = phased.Collector('Sensor', antenna, 'PropagationSpeed', c,'OperatingFrequency', fc);
receiver = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);


% Target

targetpos= [900,0,0;1000,-30,0]';
targetvel = [0,0,0;0,0,0]';

squintangle = atand(600/950);
target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', [1,1]);
pointTargets = phased.Platform('InitialPosition', targetpos,'Velocity',targetvel);

% Ground truth plot
figure(1);
h = axes;plot(targetpos(2,1),targetpos(1,1),'*b');hold on;plot(targetpos(2,2),targetpos(1,2),'*r');hold off;
set(h,'Ydir','reverse');xlim([-50 10]);ylim([800 1200]);
title('Ground Truth');ylabel('Range');xlabel('Cross-Range');

% Signal simulation using ray tracing

pm = propagationModel("raytracing", ...
    Method="sbr", ...
    CoordinateSystem="cartesian", ...
    MaxNumReflections=1, ...
    MaxNumDiffractions=0);

rxsig = zeros(truncrangesamples,numpulses);

for ii = 1:numpulses

    % Update radar platform and target position
    [radarpos, radarvel] = radarPlatform(slowTime);
    [targetpos,targetvel] = pointTargets(slowTime);

    % Define radar transmitter site
    tx = txsite("cartesian", ...
        "AntennaPosition", radarpos(:), ...
        "TransmitterFrequency", fc);

    for k = 1:size(targetpos,2)

        % Target as receiver site
        rx = rxsite("cartesian", ...
            "AntennaPosition", targetpos(:,k));

        % Ray tracing

        rays = raytrace(tx,rx,pm);

        if isempty(rays)
            continue
        end

        rayset = rays{1};
        rangeHistory(ii) = rayset(1).PropagationDistance;

        for r = 1:length(rayset)

            dist = rayset(r).PropagationDistance;
            tau = 2*dist/c;

            sig = waveform();
            sig = sig(1:truncrangesamples);
            
            delay = round(tau*fs);
            loss = db2pow(-rayset(r).PathLoss);
            phase = loss * exp(-1i*4*pi*fc*dist/c);
            
            if delay < truncrangesamples
            
                valid = 1:(truncrangesamples-delay);
            
                rxsig(delay + valid, ii) = ...
                    rxsig(delay + valid, ii) + phase*sig(valid);
            
            end
        end
    end
end

figure
plot(rangeHistory)
title('Raytrace range history')

kc = (2*pi*fc)/c;
% Compensate for the doppler due to the squint angle
rxsig=rxsig.*exp(-1i.*2*(kc)*sin(deg2rad(squintangle))*repmat(speed*eta1,1,truncrangesamples)).';

% Visualization 

imagesc(real(rxsig));title('SAR Raw Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

% Range compression

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
matchingCoeff = getMatchedFilter(waveform);
[cdata, rnggrid] = pulseCompression(rxsig, matchingCoeff);


figure;
imagesc(real(cdata));
title('SAR Range Compressed Data ');
xlabel('Cross-Range Samples');
ylabel('Range Samples');

% Azimuth compression
rma_processed = helperSquintRangeMigration(cdata,fastTime,fc,fs,prf,speed,numpulses,c,Rc,squintangle);

% Final image 

figure(2);
imagesc(abs(rma_processed(2300:3600,1100:1400).'));
title('SAR Data focused using Range Migration algorithm ') 
xlabel('Cross-Range Samples') 
ylabel('Range Samples')


%% Azimuth compression helper function 

function azcompresseddata = helperSquintRangeMigration(sigData,fastTime,fc,fs,prf,speed,numPulses,c,Rc,squintangle)
frequencyRange = linspace(fc-fs/2,fc+fs/2,length(fastTime));
krange = 2*(2*pi*frequencyRange)/c;
kaz = 2*pi*linspace(-prf/2,prf/2,numPulses)./speed;
kc = 2*pi*fc/3e8;
kazimuth = kaz.';
kus=2*(kc)*sin(deg2rad(squintangle));
kx = krange.^2-(kazimuth+kus).^2;

thetaRc = deg2rad(squintangle);
kx = sqrt(kx.*(kx > 0));

kFinal = exp(1i*(kx.*cos(thetaRc)+(kazimuth).*sin(thetaRc)).*Rc);
kfin = kx.*cos(thetaRc)+(kazimuth+kus).*sin(thetaRc);

sdata =fftshift(fft(fftshift(fft(sigData,[],1),1),[],2),2);
fsmPol = (sdata.').*kFinal;

stoltPol = fsmPol;
for i = 1:size((fsmPol),1)
    stoltPol(i,:) = interp1(kfin(i,:),fsmPol(i,:),krange(1,:));
end
stoltPol(isnan(stoltPol)) = 1e-30;

azcompresseddata = ifftshift(ifft2(stoltPol),2);

end