function pv = gpe2dpseudovorticity(gridx,gridy,psi)
    %Calculate pseudo-vorticity
    %See https://arxiv.org/pdf/1604.03595v2.pdf for details.
    dx = gridx(5)-gridx(4);
    [FXR,FYR] = gradient(real(psi),dx);
    [FXI,FYI] = gradient(imag(psi),dx);
    cp = cross([FXR(:),FYR(:),0*FXR(:)],[FXI(:),FYI(:),0*FXR(:)],2); %cross prod is 3d so imagine no gradient along z
    pv = reshape(cp(:,3),size(FXR));
end