clear all; close all; clc

% Fetching map
viewer = siteviewer(Buildings="map.osm");

% Set receiver position
rx = rxsite(Latitude=59.405220, ...
    Longitude=17.946607, ...
    AntennaHeight=30);
show(rx)

pm = propagationModel("raytracing", Method="sbr", MaxNumDiffractions=2, ...
    MaxNumReflections=3);

% Struct to store ray properties
propagationData = struct();

for i = 1:10
    
    % Moving approx. 1m per step 
    tx = txsite(Latitude=59.405023 + i*1e-5, ...
        Longitude=17.947607, ...
        TransmitterFrequency=2.5e9);
    show(tx)

    rays = raytrace(tx,rx,pm);
    
    plot(rays{1}, Type="power"); % Plotting propagation paths
    propagationData(i).Phase = NaN;
       
    if ~isempty(rays) && ~isempty(rays{1})

        % Selecting rays with specific number of interactions
        numInt = [rays{1}.NumInteractions];
        selected = rays{1}(numInt <= 3);

        %phases = [selected.PhaseShift];
        %pathLoss = [selected.PathLoss];

        % Storing phase shift and path loss for each ray
        propagationData(i).Phase = [selected.PhaseShift];
        propagationData(i).PathLoss = [selected.PathLoss];
        
        
        %idx = find([rays{1}.NumInteractions] == 1,1);

        %if ~isempty(idx)

        %propagationData(i).Phase = rays{1}(idx).PhaseShift;
        %end
    end
  
end

% Print 
for k = 1:length(rays{1})
    
    fprintf('\n========== Ray %d ==========\n', k);
    disp(rays{1}(k));
    
end

fprintf('Phase\n')
disp([propagationData.Phase])
fprintf('Path Loss\n')
disp([propagationData.PathLoss])


%% Trying to extract rays with only diffraction interactions 

selectedIdx = [];

for d = 1:length(rays{1})
    
    ray = rays{1}(d);
    
    % Extract the Type of each interaction
    interactionTypes = {ray.Interactions.Type};   % cell array of strings
    
    % Check if any of them is "Diffraction"
    if any(strcmp(interactionTypes, 'Diffraction'))
        selectedIdx(end+1) = d;
    end
    
end


diffractionRays = rays{1}(selectedIdx);

