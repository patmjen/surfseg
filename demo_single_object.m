%% Add folder with mex files
% You need to set this to your build dir
addpath surfseg/build/x64-Release/surfseg/matlab

%% Make volume containing a single ball
V = false(100,100,100);
V(50,50,50) = true;
V = single(bwdist(V) <= 30);

%% Display volume
figure(1); clf;
patch(isosurface(V,0.5),'EdgeColor','none','FaceColor',[0.8,0.8,0.8]);
axis equal;
axis([1 size(V,1) 1 size(V,1) 1 size(V,1)]);
view(45,20);
lighting phong
camlight right
material dull

%% Segment ball
% Define (very) simple region cost volume
Cost = (1 - V) - V;

% Segmentation parameters
Cen = single([50 50 50]); % Mesh center
r = single(10);           % Mesh radius
nsub = 3;                 % Subdivision level for subdivided icosahedron
nsamples = 50;            % Number of samples
step = 1;                 % Step size for samples
smoothness = 10;          % Smoothness parameter
costtype = 1;             % Cost type (1 = region costs)

[Fcs,Vtx] = mex_surfcut(Cost,Cen,r,nsub,nsamples,step,smoothness,costtype);
Fcs = Fcs{1};
Vtx = Vtx{1};
Vtx = Vtx(:,[2 1 3]) + 1; % Need to adjust vertices

%% Display segmentation
figure(2); clf;
patch(isosurface(V,0.5),'EdgeColor','none','FaceColor',[0.8,0.8,0.8],...
    'FaceAlpha',0.6);
hold on;
patch('Faces',Fcs,'Vertices',Vtx,'EdgeColor','r','FaceColor','none');

axis equal;
axis([1 size(V,2) 1 size(V,1) 1 size(V,3)]);
view(45,20);
lighting phong
camlight right
material dull
