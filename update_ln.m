function [l]=update_ln(l,th,a) %#ok<*FNDEF>
lg=logical(th<=l);
l=l.*lg+max((a*l-th)./((a-1)*l),0).*~lg;
