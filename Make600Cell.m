function [F,V] = Make600Cell(Center,rs,rt,nsubdiv)
% MAKE600CELL Make 600 cell tetrahedron mesh
% [F,V] = Make600Cell()
% [F,V] = Make600Cell(Center)
% [F,V] = Make600Cell(Center,rs)
% [F,V] = Make600Cell(Center,rs,rt)
% [F,V] = Make600Cell(Center,rs,rt,nsubdiv)
%
% Patrick M. Jensen, 2019, Technical University of Denmark

if nargin < 1
    Center = [0 0 0 0];
end
if nargin < 2
    rs = 1;
end
if nargin < 3
    rt = 1;
end
if nargin < 4
    nsubdiv = 0;
end

if iscolumn(Center)
    Center = Center';
end

V = [];
for i = 0:2^4-1
    Signs = 2*bitget(i,1:4) - 1;
    V = [V; 0.5*ones(1,4).*Signs]; %#ok<AGROW>
end
for i = 1:4
    Pi = [0 0 0 0];
    Pi(i) = 1;
    V = [V; Pi; -Pi];
end
phi = 0.5*(1 + sqrt(5));
W = 0.5*[phi,1,1/phi,0];
PermIdxs = perms(1:4);
for i = 1:size(PermIdxs,1)
    Idx = PermIdxs(i,:);
    
    % Check if permutation is even
    M = zeros(4);
    for k = 1:4
        M(k,Idx(k)) = 1;
    end
    if det(M) < 0
        % Odd permuation so skip
        continue;
    end
    
    for j = 0:2^4-1
        Signs = 2*bitget(j,1:4) - 1;
        if prod(Signs) > 0
            V = [V; W(PermIdxs(i,:)).*Signs]; %#ok<AGROW>
        end
    end
end
F = convhulln(V);

% Subdivide
for si = 1:nsubdiv
    for  i = 1:size(F,1)            
        P1 = V(F(i,1),:);
        P2 = V(F(i,2),:);
        P3 = V(F(i,3),:);
        P4 = V(F(i,4),:);
        Vn1 = 0.5*(P2 + P1);
        Vn2 = 0.5*(P3 + P1);
        Vn3 = 0.5*(P4 + P1);
        Vn4 = 0.5*(P2 + P3);
        Vn5 = 0.5*(P2 + P4);
        Vn6 = 0.5*(P3 + P4);
        V = [V; 
            Vn1/norm(Vn1); 
            Vn2/norm(Vn2); 
            Vn3/norm(Vn3);
            Vn4/norm(Vn4);
            Vn5/norm(Vn5);
            Vn6/norm(Vn6)]; %#ok<AGROW>
    end
    V = unique(V,'rows');
    F = convhulln(V);
end

% Fix orientation so all normals point 'outwards'
for i = 1:size(F,1)
    % Get tet points
    P1 = V(F(i,1),:);
    P2 = V(F(i,2),:);
    P3 = V(F(i,3),:);
    P4 = V(F(i,4),:);
    
    % Get barycenter
    PM = 0.25*(P1 + P2 + P3 + P4);
    
    % Compute edge vectors
    V1 = P2 - P1;
    V2 = P3 - P1;
    V3 = P4 - P1;
    
    % Compute normal vector, see: https://math.stackexchange.com/a/2371039
    W = [V1; V2; V3];
    N = [-det(W(:,[2 3 4])) det(W(:,[1 3 4])) -det(W(:,[1 2 4])) det(W(:,[1 2 3]))];
    
    if dot(PM,N) < 0
        % Normal points inwards so swap tet. points
        tmp = F(i,3);
        F(i,3) = F(i,4);
        F(i,4) = tmp;
    end
end

V(:,1:3) = V(:,1:3) * rs;
V(:,4) = V(:,4) * rt;
V = V + Center;