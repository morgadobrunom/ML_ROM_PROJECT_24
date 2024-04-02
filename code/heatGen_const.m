function fout = heatGen_lin(region,state,Eresults,rho,thickness, alpha)
xq = region.x;
yq = region.y;
[gradx,grady] = evaluateGradient(Eresults,xq,yq);
J_sq = (gradx./rho).^2 + (grady./rho).^2;
fout = thickness*rho.*(1 + alpha*state)*(J_sq)';