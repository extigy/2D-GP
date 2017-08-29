function [ tracks, adjacency_tracks, points] = gpeget2dvort_track(dirarg,startno,stride,endno)

n_frames = (endno-startno)/stride;
points = cell(n_frames, 1);

for i=startno:stride:endno
    [gridx,gridy,psi,pot] = gpeget2dPSI(dirarg,i);
    fprintf('read %d\n',i);
    [xlocs,ylocs,pol] = gpeget2dvort(psi,gridx,gridy,'potential',pot);
    j = (i+(stride-startno))/stride;
    points{j} = [xlocs;ylocs;pol*10]'; %*10 to make pol more "important"
    fclose('all');
end

max_linking_distance = 	6;
max_gap_closing = 24;
debug = true;
[ tracks, adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);
end