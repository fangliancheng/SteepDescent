function varargout = geodesic(z,mode)
%function of 2N variable: X1,Y1,X2,Y2...XN,YN
global numf numg 

N = 5;
a = -1;
b = -1;
c = 1; 
d = 1;
alfa = 4;
beta = 5;
delta_t = 1/(1+N);
i = (1:1:N-1)';
argout = 0;

if bitand(mode,1) 
  numf = numf + 1;
  argout = argout + 1;
  delta_t = 1/(1+N);
  ro = 1 + alfa * exp( -beta * (   z(2*i-1).^2  +  z(2*i).^2)  );
  dis = ( ( z(2*i+1) - z(2*i-1) )/delta_t ).^2 + ( ( z(2*i+2) - z(2*i) )/delta_t ).^2;
  first = ( 1 + alfa * exp( -beta*( a^2+b^2 ) ) ) * ( ( (z(1)-a)/delta_t )^2 + ( (z(2)-b)/delta_t )^2 );
  last = (1 + alfa * exp(-beta * (z(2*N-1)^2 + z(2*N)^2))) * (((c - z(2*N-1))/delta_t)^2 + ((d - z(2*N))/delta_t)^2);   
  varargout(argout) = {delta_t * ( first + sum(ro.*dis) + last)};
end

if bitand(mode,2) 
  numg = numg + 1;
  argout = argout + 1;
  DX1  = 1/delta_t*((1 + alfa * exp(-beta*(a^2+b^2)))*2*(z(1) - a) + (alfa * exp(-beta*(z(1)^2+z(2)^2))) * (-beta*2*z(1)) * ((z(3)-z(1))^2 + (z(4) - z(2))^2) + (1 + alfa * exp(-beta * (z(1)^2 + z(2)^2)))*(-2)*(z(3) - z(1)) );
  DY1  = 1/delta_t*((1 + alfa * exp(-beta*(a^2+b^2)))*2*(z(2) - b) + (alfa * exp(-beta*(z(1)^2+z(2)^2))) * (-beta*2*z(2)) * ((z(4)-z(2))^2 + (z(3)-z(1))^2) + (1 + alfa * exp(-beta * (z(1)^2 + z(2)^2)))*(-2)*(z(4) - z(2)));
  
  
  DXlast = 1/delta_t * ((1 + alfa * exp(-beta * (z(2*N-3)^2 + z(2*N-2)^2))) * 2 * (z(2*N-1) - z(2*N-3)) +...
      (alfa * exp(-beta*(z(2*N-1)^2 + z(2*N)^2))) * (-beta * 2 * z(2*N-1)) * ( (c - z(2*N-1))^2 +(d - z(2*N))^2 ) +...
      (1 + alfa * exp(-beta * (z(2*N-1)^2 + z(2*N)^2)))*(-2)*(c - z(2*N-1)));
  DYlast = 1/delta_t * ((1 + alfa * exp(-beta * (z(2*N-3)^2 + z(2*N-2)^2))) * 2 * (z(2*N) - z(2*N-2)) +...
      (alfa * exp(-beta*(z(2*N-1)^2+z(2*N)^2))) * (-beta * 2 * z(2*N)) * ( (d - z(2*N))^2 +(c - z(2*N-1))^2 ) +...
      (1 + alfa * exp(-beta * (z(2*N-1)^2 + z(2*N)^2)))*(-2)*(d - z(2*N)));

  
  for j = 2:1:N-1
      DX_mid(j) =  1/delta_t * ((1 + alfa * exp(-beta*(z(2*j-3)^2+z(2*j-2)^2)))*2*(z(2*j-1) - z(2*j-3)) + (alfa * exp(-beta*(z(2*j-1)^2+z(2*j)^2))) * (-beta*2*z(2*j-1)) * ( (z(2*j+1)-z(2*j-1))^2 + (z(2*j+2) - z(2*j))^2 )  + (1 + alfa * exp(-beta * (z(2*j-1)^2 + z(2*j)^2))) * (-2) * (z(2*j+1) - z(2*j-1)));
      DY_mid(j) =  1/delta_t * ((1 + alfa * exp(-beta*(z(2*j-3)^2+z(2*j-2)^2)))*2*(z(2*j) - z(2*j-2)) + (alfa * exp(-beta*(z(2*j-1)^2+z(2*j)^2))) * (-beta*2*z(2*j)) * ((z(2*j + 2)-z(2*j))^2 + (z(2*j+1)-z(2*j-1))^2 ) + (1 + alfa * exp(-beta * (z(2*j-1)^2 + z(2*j)^2)))*(-2)*(z(2*j+2) - z(2*j)));
  end
  gradient = [DX1;DY1];
  for k = 2:1:N-1
      gradient = [gradient;DX_mid(k);DY_mid(k)];
  end
  gradient = [gradient; DXlast; DYlast];
  varargout(argout) = {gradient};
end


return;
