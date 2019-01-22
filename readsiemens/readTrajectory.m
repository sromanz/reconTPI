function crd = readTrajectory(filename,mrprot)


ro        = mrprot.Meas.BaseResolution;
rawtrajectory = mexLoadSiemensScandata(filename);

%%
trajectory = permute(reshape(rawtrajectory,ro,3,[]),[1 3 2]);

%%
trajectory = trajectory(:,1:size(trajectory,2)/mrprot.MeasYaps.lAverages,:);

%% create arary of coordinates and normalise then in the range -0.5:0.5

crd = reshape(double(permute(trajectory,[3 1 2])),3,ro,1,[]); %% this is the data structure needed by the gridding routine

%normalise to -0.5 to 0.5 

crd(1,:)=0.5*crd(1,:)./(max(crd(1,:)));
crd(2,:)=0.5*crd(2,:)./(max(crd(2,:)));
crd(3,:)=0.5*crd(3,:)./(max(crd(3,:)));

crd(1,:)=-0.5*crd(1,:)./(min(crd(1,:)));
crd(2,:)=-0.5*crd(2,:)./(min(crd(2,:)));
crd(3,:)=-0.5*crd(3,:)./(min(crd(3,:)));

end