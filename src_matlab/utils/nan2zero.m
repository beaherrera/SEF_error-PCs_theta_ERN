function x = nan2zero(x)

x(isnan(x))=0;

end