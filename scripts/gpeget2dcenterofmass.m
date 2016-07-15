function [comx,comy] = gpeget2dcenterofmass(dens,gridx,gridy)
    dims = size(dens);
    for yy=1:(dims(2))
    for xx=1:(dims(1))
    	dens_comx(yy,xx) = dens(yy,xx)*gridx(xx);
    	dens_comy(yy,xx) = dens(yy,xx)*gridy(yy);
    end
    end
    comx = trapz(trapz(dens_comx))*(gridx(2)-gridx(1))*(gridy(2)-gridy(1));
    comy = trapz(trapz(dens_comy))*(gridx(2)-gridx(1))*(gridy(2)-gridy(1));
end