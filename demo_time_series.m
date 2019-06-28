%% Add folder with mex files
% You need to set this to your build dir
addpath surfseg/build/x64-Release/surfseg/matlab

%% Make volumetric time series with two growing balls
% This may take a few seconds...
V = ones(200,100,100,100,'single');
for t = 1:100
    V1 = false(200,100,100);
    V1(125,50,50) = true;
    V1 = single(bwdist(V1) <= 5 + round(t*(40/100)));
    
    V2 = false(200,100,100);
    V2(75,50,50) = true;
    V2 = single(bwdist(V2) <= 5 + round(t*(40/100)));
    
    V(:,:,:,t) = max(V1,V2);
end

%% Display volumetric time series
for t = 1:5:size(V,4)
    cla;    
    patch(isosurface(V(:,:,:,t),0.5),...
        'EdgeColor','none','FaceColor',[0.8,0.8,0.8]);
    axis equal;
    axis([1 size(V,2) 1 size(V,1) 1 size(V,3)]);
    title(sprintf('t = %d',t));
    view(80,10);
    lighting phong
    camlight right
    material dull
    
    pause(0.5);
end

%% Compute (very) simple region cost volume
Cost = (1 - V) - V;

%% Segment balls using 600-cell as initialization
% Segmentation parameters
Cen = [100,50,50,95]; % Mesh center in 4D space
rs = 20;              % Spatial radius
rt = 5;               % Temporal radius
nsub = 2;             % Subdivision level
nsamples = 100;       % Number of samples
step = 1;             % Step size for samples
smoothness = 40;      % Smoothness parameter
costtype = 1;         % Cost type (1 = region costs)
bend = false;         % Experimental parameter; always set this to false

[FcsI,VtxI] = Make600Cell(Cen,rs,rt,nsub); % Init mesh
[Fcs,Vtx] = mex_surfcut_4d(Cost,int32(FcsI),single(VtxI),nsamples,step,...
    smoothness,costtype,bend);

%% Display segmentation
for t = 1:5:size(V,4)
    cla;
    % Display time slice
    patch(isosurface(V(:,:,:,t),0.5),...
        'EdgeColor','none','FaceColor',[0.8,0.8,0.8],'FaceAlpha',0.6);
    hold on;
    
    % Compute and display cross section of segmentation hypersurface
    VtxT = EnsureTransversal(Vtx,Fcs,t);
    [Vs,Fs] = MeshTimeSlice(VtxT,Fcs,t);
    Vs = Vs(:,[2 1 3]) + 1; % Adjust vertices
    patch('Faces',Fs,'Vertices',Vs,'EdgeColor','r','FaceColor','none');
    
    axis equal;
    axis([1 size(V,2) 1 size(V,1) 1 size(V,3)]);
    title(sprintf('t = %d',t));
    view(80,10);
    lighting phong
    camlight right
    material dull
    
    pause(0.5);
end
