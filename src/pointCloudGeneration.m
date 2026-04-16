%% POINT CLOUD GENERATION

function points = pointCloudGeneration(x0, y0, z0, mapTarget, cubeTarget)
    points = [];

    if mapTarget
        scene = uavScenario(ReferenceLocation = [59.402979 17.956008 0]);
    
        xLim = [-100 100];
        yLim = [-100 100];
    
        addMesh(scene, "buildings", ...
            {"map2.osm", xLim, yLim, "auto"}, [0.6 0.6 0.6]);
    
        spacing = 5;   % edge sampling resolution
        radarPos = [x0, y0, z0];
        points = [];

        for k = 1:numel(scene.Meshes)
            pts = extractVisibleEdges(scene.Meshes{k}, radarPos, spacing);
            points = [points; pts];
        end
        % for k = 1:numel(scene.Meshes)
        % 
        %     verts = scene.Meshes{k}.Vertices;
        %     faces = scene.Meshes{k}.Faces;
        % 
        %     edgeSet = [];
        %     faceNormals = zeros(size(faces,1),3);
        %     faceCenters = zeros(size(faces,1),3);
        % 
        %     for f = 1:size(faces,1)
        % 
        %         v1 = verts(faces(f,1),:);
        %         v2 = verts(faces(f,2),:);
        %         v3 = verts(faces(f,3),:);
        % 
        %         n = cross(v2 - v1, v3 - v1);
        %         if norm(n) == 0
        %             continue;
        %         end
        % 
        %         n = n / norm(n);
        % 
        %         faceNormals(f,:) = n;
        %         faceCenters(f,:) = (v1 + v2 + v3) / 3;
        % 
        %         % collect triangle edges
        %         edgeSet = [edgeSet;
        %                    faces(f,1) faces(f,2);
        %                    faces(f,2) faces(f,3);
        %                    faces(f,3) faces(f,1)];
        %     end
        % 
        %     edgeSet = sort(edgeSet,2);
        %     edgeSet = unique(edgeSet,'rows');
        % 
        %     % Keep edges belonging to faces that face radar
        % 
        %     visibleEdges = [];
        % 
        %     for e = 1:size(edgeSet,1)
        % 
        %         vA = verts(edgeSet(e,1),:);
        %         vB = verts(edgeSet(e,2),:);
        % 
        %         mid = (vA + vB) / 2;
        % 
        %         r = radarPos - mid;
        % 
        %         [~, idx] = min(vecnorm(faceCenters - mid,2,2));
        % 
        %         if idx == 0 || all(faceNormals(idx,:) == 0)
        %             continue;
        %         end
        % 
        %         n = faceNormals(idx,:);
        % 
        %         if dot(n, r) > 0
        %             visibleEdges = [visibleEdges; edgeSet(e,:)];
        %         end
        %     end
        % 
        %     visibleEdges = unique(visibleEdges,'rows');
        % 
        %     for e = 1:size(visibleEdges,1)
        % 
        %         v1 = verts(visibleEdges(e,1),:);
        %         v2 = verts(visibleEdges(e,2),:);
        % 
        %         L = norm(v2 - v1);
        %         npts = max(2, ceil(L / spacing));
        % 
        %         for i = 0:npts
        %             t = i / npts;
        %             p = (1 - t)*v1 + t*v2;
        %             points = [points; p];
        %         end
        %     end
        % 
        %     if ~isempty(visibleEdges)
        %         cornerVerts = unique(visibleEdges(:));
        % 
        %         for i = 1:length(cornerVerts)
        %             v = verts(cornerVerts(i),:);
        % 
        %             for rep = 1:3
        %                 points = [points; v];
        %             end
        %         end
        %     end
        % end
    end
    
    points = unique(round(points,3),'rows');
    
    if cubeTarget
         spacing = 8;   

        cubeSize = 16;    
        halfSize = cubeSize / 2;

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

        theta = deg2rad(45);  
    
        Rz = [
            cos(theta), -sin(theta), 0;
            sin(theta),  cos(theta), 0;
            0,           0,          1
        ];
    
        verts = (Rz * verts')';
    
        edges = [
            1 2; 2 3; 3 4; 4 1;  
            5 6; 6 7; 7 8; 8 5;  
            1 5; 2 6; 3 7; 4 8   
        ];
    
        for e = 1:size(edges,1)
            v1 = verts(edges(e,1),:);
            v2 = verts(edges(e,2),:);
    
            edgeLength = norm(v2 - v1);
            npts = max(2, ceil(edgeLength / spacing));
    
            for i = 0:npts
                t = i / npts;
                p = v1 + t*(v2 - v1);
                points = [points; p];
            end
        end
    
        points = unique(points,'rows');

        % --- Center point cloud at (0,0,0) ---
        centroid = mean(points, 1);
        points = points - centroid;
    end 
      

end 