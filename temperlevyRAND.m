function randDraws = temperlevyRAND(nvals, gam, a, b)
% Range of ND sediment residence times (T)
T = logspace(log10(gam),log10(gam.*2E3),1E5); 
[F, ~] = temperlevy(T', gam, a, b); %temperLevy S(x) generation
invCDF = fit(1-F',T','linearinterp'); %inverse generation
randDraws = feval(invCDF,rand(nvals)); %random values
