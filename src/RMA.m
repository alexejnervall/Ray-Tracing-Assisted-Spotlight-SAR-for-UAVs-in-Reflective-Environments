%% Azimuth compression helper function 

function azcompresseddata = RMA(sigData,fastTime,fc,fs,prf,speed,n,c,Rc,squintAngle)

% Creating RF frequency grid for sampled signal 
frequencyRange = linspace(fc-fs/2,fc+fs/2,length(fastTime));

% Range spacial frequency axis 
krange = 2*(2*pi*frequencyRange)/c;

% Azimuth spacial frequency axis 
kaz = 2*pi*linspace(-prf/2,prf/2,n)./speed;

% Carrier spacial frequency 
kc = 2*pi*fc/3e8;
kazimuth = kaz.';

% Doppler center shift 
kus=2*(kc)*sin(deg2rad(squintAngle)); % zero for squint = 0 

% SAR geometry equation
kx = krange.^2-(kazimuth+kus).^2;

% Converting squint angle to rad
thetaRc = deg2rad(squintAngle);
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