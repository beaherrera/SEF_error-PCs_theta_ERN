function quadz = cal_quad_CSD(zs, iCSD, rc)
%CAL_QUAD_CSD calculate the current quadrupole moment from the CSD
%   calculates the current quadrupole moment from the laminar current source
%   density (CSD) relative to the center of a cylindrical cortical column
%   of 'rc' mm radius, which also corresponds to the center of 
%   the solumn's coordinate system. 
%
% Inputs:
%   zs: [mm] depth relative to the pia matter at which the CSD was
%       calculated [Nelectrodes x 1]
%   iCSD: [nA/mm3] current source density [Nelectrodes x Ntimepts]
%   z: cortical depth -> m (zero coordinate at pia matter)
%   rc: [mm] cortical column radius
% 
% Output:
%   quad: [nA*mm2] current dipole moment amplitude
%
% Author: Beatriz Herrera
%
%%

quadz = ((zs-median(zs)).^2.*mean(diff(zs)))'*iCSD*(pi*(rc^2)); % nA*mm2

end