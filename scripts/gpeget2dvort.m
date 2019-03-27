%%THE 2D VORTEX CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MIT Licence 2013-2017
%CONTRIBUTORS
%    Original implementation & current maintainer: G.Stagg
%    Bug reports & testing: P.Comaron, M.Mesgarnezhad
%    Performance improvements: E.Rickinson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xlocs,ylocs,pol] = gpeget2dvort(psi,gridx,gridy,varargin)
    dx = gridx(2)-gridx(1);
    dy = gridy(2)-gridy(1);
    dims = size(psi);
    
    %parse function arguments
    p = inputParser;
    addRequired(p,'psi');
    addRequired(p,'gridx');
    addRequired(p,'gridy');
    addParameter(p,'potential',zeros(dims));
    addParameter(p,'potentialheight',25);
    addParameter(p,'plot',0);
    addParameter(p,'plotphase',0);
    addParameter(p,'plotvortfield',0);
    addParameter(p,'boundary','periodic');
    addParameter(p,'filtersize',0);
    addParameter(p,'verbose',0);
    parse(p,psi,gridx,gridy,varargin{:});
    
    log(['Grid spacing is ',num2str(dx),',',num2str(dy),'.'], p.Results.verbose);
    log(['Box size is ',num2str(dims(2)),',',num2str(dims(1)),'.'], p.Results.verbose);
    
    %Gaussian Filter
    if p.Results.filtersize > 0
        log(['Using gaussian filter of width ',num2str(p.Results.filtersize),'.'], p.Results.verbose);
        h = fspecial('gaussian', dims, p.Results.filtersize);
        psi = imfilter(psi, h,'circular');
        phase = angle(psi);
    elseif p.Results.filtersize == 0
        log('Gaussian filter disabled.', p.Results.verbose);
        phase = angle(psi);
    else
        error('Invalid filter size.');
    end
    
    %Boundary conditions
    if strcmp(p.Results.boundary,'periodic')
        phase=[phase(end,end),phase(end,:),phase(end,1); ...
            phase(:,end),phase,phase(:,1); ...
            phase(1,end),phase(1,:),phase(1,1)];
    elseif strcmp(p.Results.boundary,'zero')
        phase=[zeros([1 dims(2)+2]); ...
            [zeros([dims(1),1]),phase,zeros([dims(1) 1])]; ...
            zeros([1 dims(2)+2])];
    elseif strcmp(p.Results.boundary,'reflective')
        phase=[phase(1,1),phase(1,:),phase(1,end); ...
            phase(:,1),phase,phase(:,end); ...
            phase(end,1),phase(end,:),phase(end,end)];
    else
        error('Invalid boundary condition.');
    end
    
    %gradient(phase) with finite differences
    velx=real(mod(phase(2:end-1,3:end)-phase(2:end-1,1:end-2)+pi,2*pi)-pi)/2;
    vely=real(mod(phase(3:end,2:end-1)-phase(1:end-2,2:end-1)+pi,2*pi)-pi)/2;
    
    %More boundary conditions
    if strcmp(p.Results.boundary,'periodic')
        velx=[velx(end,end),velx(end,:),velx(end,1); ...
            velx(:,end),velx,velx(:,1); ...
            velx(1,end),velx(1,:),velx(1,1)];
        vely=[vely(end,end),vely(end,:),vely(end,1); ...
            vely(:,end),vely,vely(:,1); ...
            vely(1,end),vely(1,:),vely(1,1)];
    elseif strcmp(p.Results.boundary,'zero')
        velx=[zeros([1 dims(2)+2]); ...
            [zeros([dims(1),1]),velx,zeros([dims(1) 1])]; ...
            zeros([1 dims(2)+2])];
        vely=[zeros([1 dims(2)+2]); ...
            [zeros([dims(1),1]),vely,zeros([dims(1) 1])]; ...
            zeros([1 dims(2)+2])];
    elseif strcmp(p.Results.boundary,'reflective')
        velx=[velx(1,1),velx(1,:),velx(1,end); ...
            velx(:,1),velx,velx(:,end); ...
            velx(end,1),velx(end,:),velx(end,end)];
        vely=[vely(1,1),vely(1,:),vely(1,end); ...
            vely(:,1),vely,vely(:,end); ...
            vely(end,1),vely(end,:),vely(end,end)];
    else
        error('Invalid boundary condition.');
    end
    
    %Get vortex field using line integrals
    l1=(velx(1:end-2,1:end-2)+velx(1:end-2,2:end-1)+velx(1:end-2,3:end));
    l2=(vely(1:end-2,1:end-2)+vely(2:end-1,1:end-2)+vely(3:end,1:end-2));
    l3=(velx(3:end,1:end-2)+velx(3:end,2:end-1)+velx(3:end,3:end));
    l4=(vely(1:end-2,3:end)+vely(2:end-1,3:end)+vely(3:end,3:end));
    presort=l4-l2+l1-l3;
    
    %Mask using density and potential
    dens = psi.*conj(psi);
    maskDens = dens>0.01*max(dens(:));
    maskDensFilled = imfill(maskDens,'holes');
    presort(~maskDensFilled)=0;
    if(sum(p.Results.potential(:)) > 0)
        presort(p.Results.potential>p.Results.potentialheight)=0;
    end
    
    if(p.Results.plotvortfield == 1)
       h=figure();
       imagesc(gridx,gridy,presort);
    end
    
    %label vort field
    negareas = bwlabel(presort<-6.2); %just under 2pi
    posareas = bwlabel(presort>6.2);

    %If bwlabel straddles a boundary pick a side and remove one of the labels.
    posareas(1,posareas(1,1:dims(2))~=0 & posareas(dims(1),1:dims(2))~=0) = 0;
    posareas(posareas(1:dims(1),1)~=0 & posareas(1:dims(1),dims(2))~=0,1) = 0;
    negareas(1,negareas(1,1:dims(2))~=0 & negareas(dims(1),1:dims(2))~=0) = 0;
    negareas(negareas(1:dims(1),1)~=0 & negareas(1:dims(1),dims(2))~=0,1) = 0;   

    %Get vortex centres
    xylocs = regionprops(posareas,'centroid');
    xylocs = cat(1,xylocs.Centroid);
    if sum(size(xylocs))>0
        xlocs = xylocs(:,1);
        ylocs = xylocs(:,2);
        pol = ones([length(xlocs),1]);
    else
        xlocs = [];
        ylocs = [];
        pol = [];
    end
    xylocs = regionprops(negareas,'centroid');
    xylocs = cat(1,xylocs.Centroid);
    if sum(size(xylocs))>0
        xlocs = [xlocs;xylocs(:,1)];
        ylocs = [ylocs;xylocs(:,2)];
        pol = [pol;repmat(-1,[length(xylocs(:,1)),1])];
    end
    %rm NaNs
    pol(isnan(xlocs))=[];
    xlocs(isnan(xlocs))=[];
    ylocs(isnan(ylocs))=[];
    %shift to grid
    xlocs = (xlocs-1)*dx + gridx(1); 
    ylocs = (ylocs-1)*dy + gridy(1);
    %delete vortices on the boundary
    if strcmp(p.Results.boundary,'zero')
        xlocs_new = xlocs(ylocs > gridy(3) & ylocs < gridy(end-3) & ylocs > gridy(3) & ylocs < gridy(end-3));
        ylocs_new = ylocs(ylocs > gridy(3) & ylocs < gridy(end-3) & ylocs > gridy(3) & ylocs < gridy(end-3));
        pol = pol(ylocs > gridy(3) & ylocs < gridy(end-3) & ylocs > gridy(3) & ylocs < gridy(end-3));
        xlocs = xlocs_new;
        ylocs = ylocs_new;
    end
    
    pos='r';
    neg='b';
    ptsize=5;
    
    if(p.Results.plot == 1)
        log('Plotting...',p.Results.verbose);
        phase=phase(2:end-1,2:end-1); %drop boundaries
        h=figure();
        if(p.Results.plotphase == 1)
            %Phase plot with brightness=density;
            C = hot;
            L = size(C,1);
            Gs = round(interp1(linspace(min(phase(:)),max(phase(:)),L),1:L,phase));
            H = reshape(C(Gs,:),[size(Gs) 3]).*repmat(mat2gray(abs(psi).^2),[1,1,3]);
            image(gridx,gridy,H);
            pos='w';
            neg='w';
            ptsize = 8;
        else
            imagesc(gridx,gridy,abs(psi).^2);
        end
        colormap(gray);
        axis image;
        axis xy;
        hold on;
        g = gscatter(xlocs,ylocs,pol,['b','r'],['^','o'],ptsize,'off');
        if(length(g)==1 && pol(1)==1)
            set(g(1), 'MarkerFaceColor', pos)
            set(g(1),'Marker','o');
            set(g(1),'MarkerEdgeColor','none');
        elseif(length(g)==1 && pol(1)==-1)
            set(g(1), 'MarkerFaceColor', neg)
            set(g(1),'Marker','^');
            set(g(1),'MarkerEdgeColor','none');
        else
            set(g(1),'MarkerEdgeColor','none');
            set(g(1), 'MarkerFaceColor', neg)
            set(g(2),'MarkerEdgeColor','none');
            set(g(2), 'MarkerFaceColor', pos)
        end
        xlabel('x', 'FontSize',16);
        ylabel('y', 'FontSize',16);
        ax = findobj(h,'type','axes','Tag','');
        set(ax,'FontSize',16)
    end
    log('Done!',p.Results.verbose);
end

function log(obj,verbose)
    if verbose > 0
        disp(obj);
    end
end