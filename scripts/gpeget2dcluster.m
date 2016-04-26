function [clusters,dipoles,free] = gpeget2dcluster(xlocs,ylocs,pol)
clusters{1} = [];
clusterso{1} = [];
dipoles{1} = [];

points =  [xlocs;ylocs]';
pold =  points;
polarity = pol';
polold = polarity;
rmdipoles = 1;
onlyoneleft=0;
founddipoles=0;
while(rmdipoles)
    %get nn for each vortex
    nnList = knnsearch(points,points,'k',2);
    if(numel(nnList)<2)
        onlyoneleft=1;
        break;
    end
    nnList = nnList(:,2);
    for i=1:length(nnList)
      if(i == nnList(nnList(i)) && polarity(i) ~= polarity(nnList(i)))
          %disp('Opposite sign mutal NN found, removing both as dipoles')
          founddipoles = founddipoles + 1;
          dipoles{founddipoles} = [points(i,:);points(nnList(i),:)];    
          points([i,nnList(i)],:) = [];
          polarity([i,nnList(i)]) = [];
          break;
      end
      if(i == length(nnList))
        %yay
        rmdipoles=0;
      end
    end
end

if(onlyoneleft==0)
    %get nn for each vortex
    nnList = knnsearch(points,points,'k',2);
    nnList = nnList(:,2);
    D = pdist2(points,points);
    for i=1:length(points)
        for j=1:length(points)
            if(i==j)
                continue;
            end;
            if(polarity(i) == polarity(j))
                 distij = D(i,j);
                 distsi = D(i,:);
                 dists_nPoli = distsi(polarity ~= polarity(i));
                 distsj = D(j,:);
                 dists_nPolj = distsj(polarity ~= polarity(j));
                 minDistTonPol = min([dists_nPoli,dists_nPolj]);

                 if(distij < minDistTonPol)
                    %disp('Clustering two vortices...')
                    clusterso=clusterup(clusterso,i,j);
                 end
            end
        end
    end
end
clusterso = clusterso(~cellfun('isempty',clusterso));  
for cl=1:length(clusterso)
    clusters{cl} = [points(clusterso{cl},:),polarity(clusterso{cl})];
end
free=findfreevorts(clusters,dipoles,pold,polold);

function [clusters]=clusterup(clusters,a,b)
    founda = 0;
    foundb = 0;
    for k=1:length(clusters)
        if (founda == 0 && sum(clusters{k} == a)>0)
            founda = k;  
        end
        if (foundb == 0 && sum(clusters{k} == b)>0)
            foundb = k;
        end
        end
        if(founda >0 && foundb == 0)
            %fprintf('Added vortex %d to cluster %d\n',b,founda);
            clusters{founda} = [clusters{founda},b];         
        end
        if(founda == 0 && foundb > 0)
            %fprintf('Added vortex %d to cluster %d\n',a,foundb);
            clusters{foundb} = [clusters{foundb},a];      
        end       
        
        if(founda > 0 && foundb > 0 && founda ~=foundb)
            %fprintf('Joining cluster %d and %d\n',founda,foundb,a,b);
            clusters{foundb} = [clusters{foundb},clusters{founda}];      
            clusters{founda} = [];
        end        
        
        if (founda == 0 && foundb == 0)
            %fprintf('Added vortex %d and %d to a new cluster, %d\n',a,b,length(clusters)+1);
            clusters{length(clusters)+1} = [a,b];
        end
    end
end

function [free]=findfreevorts(clusters,dipoles,pold,polold)
    for k=1:length(clusters)
       dims = size(clusters{k});
       for l=1:dims(1)
          for m=1:length(pold)
            if(clusters{k}(l,1) == pold(m,1) && clusters{k}(l,2) == pold(m,2))
                pold(m,:) = [];
                polold(m) = [];
                break;
            end
          end
       end
    end
    for k=1:length(dipoles)
       dims = size(dipoles{k});
       for l=1:dims(1)
          for m=1:length(pold)
            if(dipoles{k}(l,1) == pold(m,1) && dipoles{k}(l,2) == pold(m,2))
                pold(m,:) = [];
                polold(m) = [];
                break;
            end
          end
       end
    end
    free{1} = [pold,polold];
end

