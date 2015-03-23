function [drifter,dist] = gpeget2dvort_countdrift(total,tracks,adjacency_tracks, points)

n_tracks = numel(tracks);

for i = 1:length(points);
    points{i} = [points{i},repmat(i,length(points{i}),1)];
end
all_points = vertcat(points{:});

%undamped
disconsti = 9.8;
disconsto = 9.868;

drifter = zeros(length(points),1);
for i_track = 1 : n_tracks
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
 
    dists = sqrt(track_points(1,1).^2+track_points(1,2).^2);
    diste = sqrt(track_points(end,1).^2+track_points(end,2).^2);
    dist{i_track} = [dists,diste];
    
    if(dists <  disconsti && diste >  disconsti)
        if(track_points(end,4) > 85 && abs(dists-diste) > 3 )
        for i = track_points(end,4):length(points);
            drifter(i) = drifter(i) + 1;
        end
        end
    end
    if(dists >  disconsto && diste <  disconsto)
        if(track_points(end,4) > 85 && abs(dists-diste) > 3 )
        for i = track_points(end,4):length(points);
            drifter(i) = drifter(i) - 1;
        end
        end
    end
    
end


repmax = repmat(total(84),1801,1);
for i = 1:84
    drifter(i) = 0;
    repmax(i) = 0;
end
drifter(length(points))=drifter(length(points)-1)

% clf
% bar(((0:1800)*0.0005*1000)/(15*2*pi),repmax)
% hold on
% bar(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total(84)-drifter','s',15))
% bar(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total,'s',15))
% plot(((0:1800)*0.0005*1000)/(15*2*pi),repmax)
% plot(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total(84)-drifter','s',15))
% plot(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total,'s',15))
% axis([0.45, 8, 0, total(84)])

clf

drifter2 = drifter(1510)-drifter';
stuff = total(84)-drifter' - total;
stuff = max(stuff(1000:1510))-stuff;
stuff=tsmovavg(stuff(70:1510),'s',15);
drifter2=tsmovavg(drifter2(70:1510),'s',15);

total(1510)
curvefit(((84:1510)*0.0005*1000)/(15*2*pi),drifter2(15:end),stuff(15:end),total(1510));

drifter2(end)
stuff(end)
fclose('all');
end