
function [radarPosDisturbed, trajError] = trajectoryGeneration(radarPosNom, sx, sz, disturbances)

    N = size(radarPosNom, 2);
    trajError = zeros(3, N);

    dx = 0;
    dz = 0;

    if disturbances
        for ii = 1:N 

            dx = dx + sx*randn;
            dz = dz + sz*randn;
            trajError(:, ii) = [dx; 0; dz]; 
            
        end
    end 

    radarPosDisturbed = radarPosNom + trajError;

end 

   
    
