function points = extractVisibleEdges(mesh, radarPos, spacing)

    verts = mesh.Vertices;
    faces = mesh.Faces;

    nFaces = size(faces,1);
    faceNormals = zeros(nFaces,3);
    faceCenters = zeros(nFaces,3);

    for f = 1:nFaces
        v1 = verts(faces(f,1),:);
        v2 = verts(faces(f,2),:);
        v3 = verts(faces(f,3),:);

        n = cross(v2 - v1, v3 - v1);
        if norm(n) == 0
            continue;
        end

        n = n / norm(n);

        faceNormals(f,:) = n;
        faceCenters(f,:) = (v1 + v2 + v3) / 3;
    end

    edgeMap = containers.Map('KeyType','char','ValueType','any');

    for f = 1:nFaces
        vids = faces(f,:);
        edges = [vids([1 2]);
                 vids([2 3]);
                 vids([3 1])];

        for i = 1:3
            e = sort(edges(i,:));
            key = sprintf('%d_%d', e(1), e(2));

            if ~isKey(edgeMap, key)
                edgeMap(key) = [];
            end

            edgeMap(key) = [edgeMap(key), f];
        end
    end

    visibleEdges = [];

    keysList = keys(edgeMap);

    for i = 1:length(keysList)

        key = keysList{i};
        faceIdx = edgeMap(key);
        vertsIdx = sscanf(key, '%d_%d')';

        keepEdge = false;

        if length(faceIdx) == 1
            keepEdge = true;

        else
            n1 = faceNormals(faceIdx(1),:);
            n2 = faceNormals(faceIdx(2),:);

            if norm(n1) == 0 || norm(n2) == 0
                continue;
            end

            if dot(n1, n2) > 0.99
                continue;
            end

            keepEdge = true;
        end

        if ~keepEdge
            continue;
        end

        vA = verts(vertsIdx(1),:);
        vB = verts(vertsIdx(2),:);
        mid = (vA + vB) / 2;

        r = radarPos - mid;

        isVisible = false;

        for f = faceIdx
            n = faceNormals(f,:);
            if norm(n) == 0
                continue;
            end

            if dot(n, r) > 0
                isVisible = true;
                break;
            end
        end

        if isVisible
            visibleEdges = [visibleEdges; vertsIdx];
        end
    end

    visibleEdges = unique(visibleEdges,'rows');

    points = [];

    for e = 1:size(visibleEdges,1)

        v1 = verts(visibleEdges(e,1),:);
        v2 = verts(visibleEdges(e,2),:);

        L = norm(v2 - v1);
        npts = max(2, ceil(L / spacing));

        for i = 0:npts
            t = i / npts;
            p = (1 - t)*v1 + t*v2;
            points = [points; p];
        end
    end

    if ~isempty(visibleEdges)
        cornerVerts = unique(visibleEdges(:));

        for i = 1:length(cornerVerts)
            v = verts(cornerVerts(i),:);

            for rep = 1:3
                points = [points; v];
            end
        end
    end
end