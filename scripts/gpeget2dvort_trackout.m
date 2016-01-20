function [ tracks, adjacency_tracks, points] = gpeget2dvort_trackout(dirarg,startno,stride,endno,speed,nx,ny)

n_frames = (endno-startno)/stride;
points = cell(n_frames, 1);

for i=startno:stride:endno
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    [comx,comy] = gpegetcenterofmass(dens,gridx,gridy);
    [xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);
    j = (i+(stride-startno))/stride;
    points{j} = [xlocs-comx;ylocs-comy;pol*10]';
end

max_linking_distance = 	2;
max_gap_closing = 4;
debug = true;
[ tracks, adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);


fclose('all');
end