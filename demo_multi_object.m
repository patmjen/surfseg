%% Add folder with mex files
% You need to set this to your build dir
addpath surfseg/build/x64-Release/surfseg/matlab

%% Make volume containing three overlapping balls
V1 = false(100,100,100);
V1(35,50,35) = true;
V1 = single(bwdist(V1) <= 30);

V2 = false(100,100,100);
V2(70,50,35) = true;
V2 = single(bwdist(V2) < 30);

V3 = false(100,100,100);
V3(50,50,70) = true;
V3 = single(bwdist(V3) < 30);

V = max(V1,V2);
V = max(V,V3);
clear V1 V2 V3;

%% Display volume
figure(1); clf;
patch(isosurface(V,0.5),'EdgeColor','none','FaceColor',[0.8,0.8,0.8]);
axis equal;
axis([1 size(V,2) 1 size(V,1) 1 size(V,3)]);
view(45,20);
lighting phong
camlight right
material dull

%% Segment balls
% Define (very) simple region cost volume
Cost = (1 - V) - V;

% Segmentation parameters
Cen = single([35 50 35;             % Mesh centers
              70 50 35;
              50 50 70]);
R = single([10,10,10]);             % Mesh radii
nsub = 3;                           % Subdiv. lvl. for subdiv. icosahedron
nsamples = 50;                      % Number of samples
step = 1;                           % Step size for samples
smoothness = 10;                    % Smoothness parameter
costtype = 1;                       % Cost type (1 = region costs)
Conn = ones(3) - eye(3);            % Mesh connectivity
Conn = logical(sparse(Conn));

[Fcs,Vtx] = mex_surfcut_planesep_qpbo(Cost,Cen,R,nsub,nsamples,step,...
    smoothness,Conn,costtype);

%% Display segmentation
figure(2); clf;
patch(isosurface(V,0.5),'EdgeColor','none','FaceColor',[0.8,0.8,0.8],...
    'FaceAlpha',0.6);
hold on;
Colors = {'r','g','b'};
for i = 1:length(Fcs)
    % Note that vertices need to be adjusted
    patch('Faces',Fcs{i},'Vertices',Vtx{i}(:,[2 1 3])+1,...
        'EdgeColor',Colors{i},'FaceColor','none');
end

axis equal;
axis([1 size(V,1) 1 size(V,1) 1 size(V,1)]);
view(45,20);
lighting phong
camlight right
material dull
