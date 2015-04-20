function [FX0,FX1,FX2,FX3,FY0,FY1,FY2,FY3] = gpeget2dforcesplit(psi,gridx,gridy,potential)
    dx = gridx(5)-gridx(4);
    dims = size(real(psi));
   
    density = psi.*conj(psi);
    phase = atan2(imag(psi),real(psi));
    [~,~,velocity(1,:,:),velocity(2,:,:)] = gpe2dvelocity(gridx,gridy,phase);
    for j=1:2,
    for k=1:2,
        T0(j,k,:,:) = density.*squeeze(velocity(j,:,:)).*squeeze(velocity(k,:,:)); 
    end
    end
    %imagesc(gridx,gridy,squeeze(T0(1,2,:,:)))
    
    pressure = density.^2;
    for j=1:2,
    for k=1:2,
        T1(j,k,:,:) = pressure.*(j==k); 
    end
    end
        
    [Dlndens(1,:,:),Dlndens(2,:,:)] = gradient(log(density+1e-14),dx);
    [D2lndens(1,1,:,:),D2lndens(2,1,:,:)] = gradient(squeeze(Dlndens(1,:,:)),dx);
    [D2lndens(1,2,:,:),D2lndens(2,2,:,:)] = gradient(squeeze(Dlndens(2,:,:)),dx);
    
    for j=1:2,
    for k=1:2,
        T2(j,k,:,:) = -(1.0/4.0).*density.*squeeze(D2lndens(j,k,:,:));
    end
    end
    
    [DT0(1,1,1,:,:),~] = gradient(squeeze(T0(1,1,:,:)),dx);
    [DT0(1,1,2,:,:),~] = gradient(squeeze(T0(1,2,:,:)),dx);
    [~,DT0(2,2,1,:,:)] = gradient(squeeze(T0(2,1,:,:)),dx);
    [~,DT0(2,2,2,:,:)] = gradient(squeeze(T0(2,2,:,:)),dx);
    [DT1(1,1,1,:,:),~] = gradient(squeeze(T1(1,1,:,:)),dx);
    [DT1(1,1,2,:,:),~] = gradient(squeeze(T1(1,2,:,:)),dx);
    [~,DT1(2,2,1,:,:)] = gradient(squeeze(T1(2,1,:,:)),dx);
    [~,DT1(2,2,2,:,:)] = gradient(squeeze(T1(2,2,:,:)),dx);
    [DT2(1,1,1,:,:),~] = gradient(squeeze(T2(1,1,:,:)),dx);
    [DT2(1,1,2,:,:),~] = gradient(squeeze(T2(1,2,:,:)),dx);
    [~,DT2(2,2,1,:,:)] = gradient(squeeze(T2(2,1,:,:)),dx);
    [~,DT2(2,2,2,:,:)] = gradient(squeeze(T2(2,2,:,:)),dx);     

    for k=1:2,
        DsumT0(k,:,:) = squeeze(DT0(1,1,k,:,:)+DT0(2,2,k,:,:));
        DsumT1(k,:,:) = squeeze(DT1(1,1,k,:,:)+DT1(2,2,k,:,:));
        DsumT2(k,:,:) = squeeze(DT2(1,1,k,:,:)+DT2(2,2,k,:,:));
    end

    [DPot(1,:,:),DPot(2,:,:)] = gradient(potential,dx);
    for k=1:2,
        rhoDv(k,:,:) = density.*squeeze(DPot(k,:,:));
    end 
    
    %imagesc(gridx,gridy,squeeze(DsumT0(1,:,:)))
    FX0 = dx*dx*simpson(simpson(squeeze(DsumT0(1,:,:)).*(potential>1)));
    FX1 = dx*dx*simpson(simpson(squeeze(DsumT1(1,:,:)).*(potential>1)));
    FX2 = dx*dx*simpson(simpson(squeeze(DsumT2(1,:,:)).*(potential>1)));
    FX3 = dx*dx*simpson(simpson(squeeze(rhoDv(1,:,:))));
    FY0 = dx*dx*simpson(simpson(squeeze(DsumT0(2,:,:)).*(potential>1)));
    FY1 = dx*dx*simpson(simpson(squeeze(DsumT1(2,:,:)).*(potential>1)));
    FY2 = dx*dx*simpson(simpson(squeeze(DsumT2(2,:,:)).*(potential>1)));
    FY3 = dx*dx*simpson(simpson(squeeze(rhoDv(2,:,:))));
    
end