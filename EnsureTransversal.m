function V = EnsureTransversal(V,T,t,tolfac)
% ENSURETRANSVERSAL Ensure mesh is transversal with given time plane
% V = EnsureTransversal(V,T,t)
% V = EnsureTransversal(V,T,t,tolfac)
%
% INPUT
% V      - N x 4 array with vertex coordinates
% T      - M x 4 array with tet. indices
% t      - Time coordinate for time plane
% tolfac - (Optional, default: 10) Tolerance factor.
%
% OUTPUT
% V      - N x 4 array with corrected vertex coordinates
%
% See also: MESHTIMESLICE
% 
% Patrick M. Jensen, 2019, Technical University of Denmark

if nargin < 4, tolfac = 10; end
tol = tolfac*eps(class(V));

for ti = 1:size(T,1)
    TIdx = T(ti,:);
    
    Signs = RobustSign(V(TIdx,4) - t, tol);
    nzer = nnz(Signs == 0);
    
    % Check for problematic intersections
    if nzer == 4
        % Intersection is entire tet
        V(TIdx,4) = V(TIdx,4) + 2*tolfac*eps(V(TIdx(1),4));
    elseif nzer == 2
        % Intersection is a line segment or triangle. If triangle, then the
        % intersecting edge might be part of a non-manifold region, so we
        % perturb the endpoints out anyway
        Idx = find(Signs == 0);
        offset = 2*tolfac*eps(V(TIdx(Idx),4));
        V(TIdx(Idx),4) = V(TIdx(Idx),4) + offset;
    elseif nzer == 1 
        % Intersection is a point or a triangle. If triangle, then 
        % neighborhood might have 'negative' curvature, so we perturb the
        % vertex out anyway
        Idx = find(Signs == 0);
        offset = 2*tolfac*eps(V(TIdx(Idx),4));
        V(TIdx(Idx),4) = V(TIdx(Idx),4) + offset;
    end
end