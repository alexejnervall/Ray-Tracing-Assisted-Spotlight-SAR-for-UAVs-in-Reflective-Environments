%% POINT CLOUD GENERATION

function points = pointCloudGeneration(x0, y0, z0)


    scene = uavScenario(ReferenceLocation = [59.402979 17.956008 0]);
    
    xLim = [-100 100];
    yLim = [-100 100];
    
    addMesh(scene,"buildings",{"map2.osm",xLim,yLim,"auto"},[0.6 0.6 0.6]); 
    
    points = [];
    
    spacing = 6;   % controls density on surfaces
    
    radarPos = [x0, y0, z0];
    
    for k = 1:numel(scene.Meshes)
    
        verts = scene.Meshes{k}.Vertices;
        faces = scene.Meshes{k}.Faces;
    
        for f = 1:size(faces,1)
    
            v1 = verts(faces(f,1),:);
            v2 = verts(faces(f,2),:);
            v3 = verts(faces(f,3),:);
    
            % --- Compute face normal ---
            normal = cross(v2 - v1, v3 - v1);
            if norm(normal) == 0
                continue;
            end
            normal = normal / norm(normal);
    
            % --- Face center ---
            center = (v1 + v2 + v3)/3;
    
            % --- Vector toward radar ---
            r = radarPos - center;
    
            % --- Keep ALL faces facing radar (any angle) ---
            if dot(normal, r) <= 0
                continue; % backface → skip
            end
    
            % --- Sample triangle surface ---
            % Estimate number of samples based on triangle size
            area = 0.5 * norm(cross(v2 - v1, v3 - v1));
            npts = max(3, ceil(area / (spacing^2)));
    
            for i = 1:npts
                % Random barycentric sampling (uniform over triangle)
                a = rand;
                b = rand;
    
                if a + b > 1
                    a = 1 - a;
                    b = 1 - b;
                end
    
                p = v1 + a*(v2 - v1) + b*(v3 - v1);
                points = [points; p];
            end
    
        end
    end
    
    % Remove duplicates
    points = unique(points,'rows');
    

    % 
    % points = [];
    % 
    % spacing = 8;   % controls point density along edges
    % 
    % cubeSize = 16;    
    % halfSize = cubeSize / 2;
    % 
    % % --- Define cube vertices ---
    % verts = [
    %     -halfSize, -halfSize, -halfSize;
    %      halfSize, -halfSize, -halfSize;
    %      halfSize,  halfSize, -halfSize;
    %     -halfSize,  halfSize, -halfSize;
    %     -halfSize, -halfSize,  halfSize;
    %      halfSize, -halfSize,  halfSize;
    %      halfSize,  halfSize,  halfSize;
    %     -halfSize,  halfSize,  halfSize
    % ];
    % 
    % theta = deg2rad(45);  % 45 degrees
    % 
    % Rz = [
    %     cos(theta), -sin(theta), 0;
    %     sin(theta),  cos(theta), 0;
    %     0,           0,          1
    % ];
    % 
    % % Rotate all vertices
    % verts = (Rz * verts')';
    % 
    % % --- Define cube edges (pairs of vertex indices) ---
    % edges = [
    %     1 2; 2 3; 3 4; 4 1;  % bottom square
    %     5 6; 6 7; 7 8; 8 5;  % top square
    %     1 5; 2 6; 3 7; 4 8   % vertical edges
    % ];
    % 
    % for e = 1:size(edges,1)
    %     v1 = verts(edges(e,1),:);
    %     v2 = verts(edges(e,2),:);
    % 
    %     % --- Compute number of points along edge ---
    %     edgeLength = norm(v2 - v1);
    %     npts = max(2, ceil(edgeLength / spacing));
    % 
    %     for i = 0:npts
    %         t = i / npts;
    %         p = v1 + t*(v2 - v1);
    %         points = [points; p];
    %     end
    % end
    % 
    % % Remove duplicates
    % points = unique(points,'rows');
    % --- Center point cloud at (0,0,0) ---


    centroid = mean(points, 1);
    points = points - centroid;

end 