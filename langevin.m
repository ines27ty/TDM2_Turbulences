function X = Langevin(Xmean, Xvar, T, dt, N)
	
	%return a signal given by the Langevin process
	% with:
	% * Xmean: the mean of the process
	% * Xvar: its variance
	% * T:its correlation time 
	% * dt: the time step
	% and N the number of time step
	
	dt_adim=dt/T;
	h=sqrt(Xvar*dt_adim);
	
	X=zeros(N,1);
	X(1)=randn()*sqrt(Xvar);
	for i=2:N
		dx = -(X(i-1) - Xmean) * dt_adim;
		dx = dx + randn()* h;
		X(i) = X(i-1) + dx ;
	end;
	
end
	 