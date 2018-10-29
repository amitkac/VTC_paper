function y_out=OutageProb(Noise_var,num,lambda_11,sl,psi,z)

lam0=lambda_11;
  lam1=sl;
  alpha=num;
  N=Noise_var;
  q=psi;
  i=1;
  
 for  zz=z
  
 fxn1= @ (zz) exp(-N.*lam1.*zz).* ((lam0+lam1.*zz).^(-alpha))...
       .* (N.*igamma(alpha,q.*(lam0+lam1.*zz))+((lam0+lam1.*zz).^(-1))...
       .*igamma(alpha+1,q.*(lam0+lam1.*zz)));
  I0=integral(fxn1,0,zz);
 
   
a=1- ((lam0.^alpha).* ((lam0+lam1.*zz).^(-alpha)).*exp(-N.*lam1.*zz));

c=lam1.*lam0.^alpha.*I0./igamma(alpha,0);

d=(1- exp(-lam1.*(N+q).*zz)).*igamma(alpha,q.*lam0)./igamma(alpha,0);

O2(i)=a-c+d;
  i=i+1;
 end
 
y_out=O2;

end