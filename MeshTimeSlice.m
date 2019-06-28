function [Vs,Fs,TIs,VNs] = MeshTimeSlice(V,T,t,trimesh,VN)
% MESHTIMESLICE Find 3D cross section of 4D mesh at a given time
% [Vs,Fs,TIs,VNs] = MeshTimeSlice(V,T,t)
% [Vs,Fs,TIs,VNs] = MeshTimeSlice(V,T,t,trimesh)
% [Vs,Fs,TIs,VNs] = MeshTimeSlice(V,T,t,trimesh,VN)
%
% INPUT
% V           - N x 4 array with vertex coodinates
% T           - M x 4 array with tet. indices
% t           - Time coordinate of time slice
% trimesh     - (Optional, default: false) Ensure output is a triangle mesh
% VN          - (Optional) N x 4 array with vertex normals. Only needed if
%               VNs output is assigned.
%
% OUTPUT 
% Vs  - Ns x 3 array with vertex corrdinates for intersection
% Fs  - Ms x (3 or 4) array with polygon indices for intersection
% TIs - Array with indices of intersecting tetrahedra
% VNs - Ns x 3 array with projected vertex normals.
%
% Patrick M. Jensen, 2019, Technical University of Denmark

if nargin < 4
    trimesh = false;
end

compNorms = nargout > 3;

nvert = size(V,1);
ntet = size(T,1);

PEMap = zeros(nvert,nvert);

try
    tol = 10*eps(class(V));
catch 
    tol = 0;
end

Vs = [];
Fs = [];
TIs = [];

if compNorms
    VNs = [];
end
for ti = 1:ntet
    % Check for intersection with time slice
    P = V(T(ti,:),:);
    
    Signs = RobustSign(P(:,4) - t, tol);
    npos = nnz(Signs > 0);
    nneg = nnz(Signs < 0);
    nzer = nnz(Signs == 0);
    if nzer == 4
        % Whole tet is in time slice so add all points        
        i1 = AddPt(T(ti,1));
        i2 = AddPt(T(ti,2));
        i3 = AddPt(T(ti,3));
        i4 = AddPt(T(ti,4));
        TIs = [TIs ti]; %#ok<AGROW>
        if trimesh
            Fs = [Fs;
                i1 i2 i3;
                i1 i2 i4;
                i1 i3 i4;
                i2 i3 i4]; %#ok<AGROW>
        else
            Fs = [Fs;
                i1 i2 i3 -1;
                i1 i2 i4 -1;
                i1 i3 i4 -1;
                i2 i3 i4 -1]; %#ok<AGROW>
        end
    elseif nzer == 0 && npos == 2 && nneg == 2
        % Intersection is a quad
        IdxPos = find(Signs > 0);
        IdxNeg = find(Signs < 0);
        [P1,i1] = GetInt(T(ti,IdxPos(1)),T(ti,IdxNeg(1)));
        [P2,i2] = GetInt(T(ti,IdxPos(1)),T(ti,IdxNeg(2)));
        [P3,i3] = GetInt(T(ti,IdxPos(2)),T(ti,IdxNeg(1)));
        [P4,i4] = GetInt(T(ti,IdxPos(2)),T(ti,IdxNeg(2)));
        Cen = 0.25*(P1 + P2 + P3 + P4);
        N = cross(P1 - Cen,P2 - Cen);
        
        % Manual sort
        if CWLess(P2,P1,Cen,N)
            tmp = i1;
            i1 = i2;
            i2 = tmp;
        end
        if CWLess(P3,P4,Cen,N)
            tmp = i3;
            i3 = i4;
            i4 = tmp;
        end
        
        TIs = [TIs ti]; %#ok<AGROW>
        if trimesh
            Fs = [Fs; i1,i2,i3; i2,i3,i4]; %#ok<AGROW>
        else
            Fs = [Fs; i1,i2,i4,i3]; %#ok<AGROW>
        end
    elseif nzer == 3
        % Triangle face intersects slice
        Idx = 1:4;
        Idx(Signs ~= 0) = [];
        
        i1 = AddPt(T(ti,Idx(1)));
        i2 = AddPt(T(ti,Idx(2)));
        i3 = AddPt(T(ti,Idx(3)));
        TIs = [TIs ti]; %#ok<AGROW>
        if trimesh
            Fs = [Fs; i1 i2 i3]; %#ok<AGROW>
        else
            Fs = [Fs; i1 i2 i3 -1]; %#ok<AGROW>
        end
    elseif npos == 1 || nneg == 1
        % Intersection is a triangle
        if npos == 1
            IdxA = find(Signs > 0);
            IdxB = find(Signs <= 0);
        else
            IdxB = find(Signs >= 0);
            IdxA = find(Signs < 0);
        end
        [~,i1] = GetInt(T(ti,IdxA(1)),T(ti,IdxB(1)));
        [~,i2] = GetInt(T(ti,IdxA(1)),T(ti,IdxB(2)));
        [~,i3] = GetInt(T(ti,IdxA(1)),T(ti,IdxB(3)));
        
        TIs = [TIs ti]; %#ok<AGROW>
        if trimesh
            Fs = [Fs; i1 i2 i3]; %#ok<AGROW>
        else
            Fs = [Fs; i1 i2 i3 -1]; %#ok<AGROW>
        end
    elseif nzer == 2 && (npos == 2 || nneg == 2)
        % Intersection is a line segment
        Idx = find(Signs == 0);
        i1 = AddPt(T(ti,Idx(1)));
        i2 = AddPt(T(ti,Idx(2)));
        TIs = [TIs ti]; %#ok<AGROW>
        if trimesh
            Fs = [Fs; i1 i2 -1]; %#ok<AGROW>
        else
            Fs = [Fs; i1 i2 -1 -1]; %#ok<AGROW>
        end
    elseif nzer == 1 && (npos == 3 || nneg == 3)
        % Intersection is a point
        TIs = [TIs ti]; %#ok<AGROW>
        Idx = find(Signs == 0);
        AddPt(T(ti,Idx(1)));
        % No edges to add
    end
    % Else there is no intersection
end
% Remove duplicate faces
[~,FIdx] = unique(sort(Fs,2),'rows');
Fs = Fs(FIdx,:);
Fs(Fs == -1) = nan;

    function i = AddPt(vi)
        if PEMap(vi,vi) == 0
            Vs = [Vs; V(vi,1:3)];
            if compNorms
                VNs = [VNs; VN(vi,:)];
            end
            i = size(Vs,1);
            PEMap(vi,vi) = i;
        else
            i = PEMap(vi,vi);
        end 
    end

    function [P,i] = GetInt(i1,i2)
        if V(i1,4) == t
            if PEMap(i1,i1) == 0
                P = V(i1,1:3);
                Vs = [Vs; P];
                if compNorms
                    VNs = [VNs; VN(i1,:)];
                end
                i = size(Vs,1);
                PEMap(i1,i1) = i;
            else
                i = PEMap(i1,i1);
                P = Vs(i,:);
            end
        elseif V(i2,4) == t
            if PEMap(i2,i2) == 0
                P = V(i2,1:3);
                Vs = [Vs; P];
                if compNorms
                    VNs = [VNs; VN(i2,:)];
                end
                i = size(Vs,1);
                PEMap(i2,i2) = i;
            else
                i = PEMap(i2,i2);
                P = Vs(i,:);
            end
        else
            if PEMap(i1,i2) == 0
                a = (V(i1,4) - t)/(V(i1,4) - V(i2,4));
                P = (1 - a)*V(i1,1:3) + a*V(i2,1:3);
                Vs = [Vs; P];
                if compNorms
                    N =  (1 - a)*VN(i1,:) + a*VN(i2,:);
                    VNs = [VNs; N];
                end
                i = size(Vs,1);
                PEMap(i1,i2) = i;
                PEMap(i2,i1) = i;
            else
                i = PEMap(i1,i2);
                P = Vs(i,:);
            end
        end
    end
end

function l = CWLess(P1,P2,Cen,N)
l = dot(N, cross(P1 - Cen, P2 - Cen)) > 0;
end