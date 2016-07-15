function track_points_all = gpe2dplottracks(tracks,adjacency_tracks, points)

n_tracks = numel(tracks);

all_points = vertcat(points{:});
track_points_all = cell(2,1);
ci = 1;
cc=jet(4);
for i_track = 1 : n_tracks
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
    if(length(track_points(:,1)) > 1250)
        plot(track_points(:,1)',track_points(:,2)','LineWidth',2,'color',cc(ci,:));
        track_points_all{1,ci} = track_points(:,1)';
        track_points_all{2,ci} = track_points(:,2)';
        ci = ci + 1;
    end
    hold on
end
end
