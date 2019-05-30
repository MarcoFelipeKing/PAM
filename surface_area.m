function A=surface_area(r,c)

m=7; v=1.9^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
A=lognrnd(mu,sigma,r,c);
end
