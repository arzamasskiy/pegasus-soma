fig(1);clf;

npts = 100;
Tratio = linspace(0.1,2.5,npts);
x = zeros(1,npts);

xo = complex(sqrt(1/(2*Tratio(1))),0);
tol  = 1e-6;

for j=1:npts
    
  tau = Tratio(j);
  
  notdone = true;
  while notdone
    
    Z = 1i*sqrt(pi)*exp(-xo^2)*(1+erfz(1i*xo));
    dZdx = -2*(1+xo*Z);
    
    f = tau - 0.5*dZdx;
    dfdx = Z + xo*dZdx;
    
    xn = xo - f/dfdx;
    
    if abs(xn-xo) > tol
      xo = xn;
    else
      xo = xn;
      w(j) =  real(xn);
      x(j) = -imag(xn);
      notdone = false;
    end
    
  end

end

conv = pi/8*sqrt(2);

semilogy(Tratio,x);
axis([0.1 2.5 1e-2 1]);
hold on;
semilogy([1/3],[0.156]/conv,'ok');
hold off;