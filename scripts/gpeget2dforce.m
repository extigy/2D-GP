function [FX, FY] = gpeget2dforce(psi,gridx,gridy,potential)
    dx = gridx(5)-gridx(4);
    dims = size(real(psi));

    [DPsi(1,:,:),DPsi(2,:,:)] = gradient(psi,dx);
    [D2Psi(1,1,:,:),D2Psi(2,1,:,:)] = gradient(squeeze(DPsi(1,:,:)),dx);
    [D2Psi(1,2,:,:),D2Psi(2,2,:,:)] = gradient(squeeze(DPsi(2,:,:)),dx);
    
    
    for j=1:2,
    for k=1:2,
        T(j,k,:,:) = real((1.0/4.0)*(conj(squeeze(DPsi(j,:,:))).*squeeze(DPsi(k,:,:))  ...
                                 -conj(psi(:,:)).*squeeze(D2Psi(j,k,:,:))         ...
                                 +squeeze(DPsi(j,:,:)).*conj(squeeze(DPsi(k,:,:)))...
                                 -psi(:,:).*conj(squeeze(D2Psi(j,k,:,:))))        ...
                    +((j==k)/2)*(psi.*conj(psi).*psi.*conj(psi)));
    end
    end

    %imagesc(gridx,gridy,real(squeeze(T(1,2,:,:))))


    [DT(1,1,1,:,:),~] = gradient(squeeze(T(1,1,:,:)),dx);
    [DT(1,1,2,:,:),~] = gradient(squeeze(T(1,2,:,:)),dx);
    [~,DT(2,2,1,:,:)] = gradient(squeeze(T(2,1,:,:)),dx);
    [~,DT(2,2,2,:,:)] = gradient(squeeze(T(2,2,:,:)),dx);

    for k=1:2,
        DsumT(k,:,:) = squeeze(DT(1,1,k,:,:)+DT(2,2,k,:,:));
    end

    %imagesc(gridx,gridy,squeeze(DsumT(2,:,:)));
    [DPot(1,:,:),DPot(2,:,:)] = gradient(potential,dx);
    
    for k=1:2,
        rhoDv(k,:,:) = psi.*conj(psi).*squeeze(DPot(k,:,:));
    end 
    %imagesc(gridx,gridy,squeeze(DsumT(2,:,:)).*(potential>1)+squeeze(rhoDv(2,:,:)));
    
    FX = dx*dx*simpson(simpson(squeeze(DsumT(1,:,:)).*(potential>1)))+dx*dx*simpson(simpson(squeeze(rhoDv(1,:,:))));
    FY = dx*dx*simpson(simpson(squeeze(DsumT(2,:,:)).*(potential>1)))+dx*dx*simpson(simpson(squeeze(rhoDv(2,:,:))));
    
end