function fr = firingrate(psth, w, flt)
% FIRINGRATE: generate continous firing rate by convolving the psth with a
% Gaussian of 'w' width.
%
% Input:
%       psth: peri-stimulus time histogram
%       w: width of the Gaussian filter
% Output:
%       fr: continous firing rate
%
%%

if flt 
    fr = smoothdata(psth,'gaussian',w);
else
    Gauss_width = 101;
    kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,w);
    fr         = conv(psth,kernel,'same');
end

end