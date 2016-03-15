function gpe2dplotclusters(clusters,dipoles,free)

%plot free vortices
plotvorts(free{1}(:,1),free{1}(:,2),free{1}(:,3));
hold on;

%plot clusters
for cnum = 1:length(clusters);
    cluster = clusters{cnum};
    dims = size(cluster);
    numVerts = dims(1);

    %make adjecency matrix
    X = ones(numVerts);
    X(logical(eye(size(X)))) = 0;

    %weights matrix
    w = pdist2(cluster(:,1:2),cluster(:,1:2));

    [~, ST, X_st] = kruskal(X, w);
    STdims = size(ST);
    %draw cluster vorts
    plotvorts(cluster(:,1),cluster(:,2),cluster(:,3));

    %draw lines
    if(cluster(1,3) == 1)
        col = 'r';
    else 
        col = 'b';
    end
    for linenum = 1:STdims(1);
        plot([cluster(ST(linenum,1),1),cluster(ST(linenum,2),1)],[cluster(ST(linenum,1),2),cluster(ST(linenum,2),2)],col);
    end
end   

%plot dipoles
for dnum = 1:length(dipoles);
    dipole = dipoles{dnum};
    %draw dipole vorts
    g = plotvorts(dipole(:,1),dipole(:,2),[1,1]);
    set(g(1), 'MarkerFaceColor', 'g')
    set(g(1),'Marker','o');
    set(g(1),'MarkerEdgeColor','none');
    %draw line
    plot([dipole(1,1),dipole(2,1)],[dipole(1,2),dipole(2,2)],'g');
end  
end

function g = plotvorts(xlocs,ylocs,pol)
    g = gscatter(xlocs,ylocs,pol,['b','r'],['^','o'],5,'off');
    if(length(g)==1 && pol(1) == 1)
        set(g(1), 'MarkerFaceColor', 'r')
        set(g(1),'Marker','o');
        set(g(1),'MarkerEdgeColor','none');
    end
    if(length(g)==1 && pol(1) == -1)
        set(g(1), 'MarkerFaceColor', 'b')
        set(g(1),'Marker','^');
        set(g(1),'MarkerEdgeColor','none');
    end
    if(length(g)==2)
        set(g(1),'MarkerEdgeColor','none');
        set(g(1), 'MarkerFaceColor', 'b')
        set(g(2),'MarkerEdgeColor','none');
        set(g(2), 'MarkerFaceColor', 'r')
    end
end