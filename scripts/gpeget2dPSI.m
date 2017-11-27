function [gridx,gridy,psi,potential] = gpeget2dPSI(dirarg,frame,varargin)
p = inputParser;
addRequired(p,'dirarg');
addRequired(p,'frame');
addParameter(p,'prefix','psi');
parse(p,dirarg,frame,varargin{:});
dirarg = regexprep(dirarg, '/$', '');
datalocation = strcat(dirarg, '/',p.Results.prefix,'.%06d.nc');
fname = sprintf(datalocation,frame);
gridx = ncread(fname,'x');
gridy = ncread(fname,'y');
real = ncread(fname,'real');
imag = ncread(fname,'imag');
potential = ncread(fname,'pot');
psi = real + 1i.*imag;
potential = potential';
psi = psi';
fclose('all');
end