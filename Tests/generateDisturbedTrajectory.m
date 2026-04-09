
function [radarPosDisturbed, trajError] = generateDisturbedTrajectory(radarPosNominal,sigma_x,sigma_z,mode)

%mode: if true Brownian noise
N = size(radarPosNominal,2);
trajError = zeros(3, N); % Initialize trajectory error matrix

dx = 0;
dz = 0;

for ii = 1:N

    if mode
        dx = dx + sigma_x*randn;
        dz = dz + sigma_z*randn;
    else
        dx = 0; % No disturbance in x
        dz = 0; % No disturbance in z
    end
    trajError(:,ii)=[dx;0;dz]; % Calculate trajectory error
    
end

radarPosDisturbed = radarPosNominal + trajError;


end