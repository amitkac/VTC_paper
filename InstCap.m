function y=InstCap(Noise_var,num,lambda_11,sl,psi,z)

N=Noise_var;
alpha=num;
beta=lambda_11;
z1=exp(z)-1;
ss= alpha+ N.*(beta+sl.*z1).*(1-igamma(alpha,psi.*(beta+sl.*z1))...
       ./igamma(alpha,0))-(igamma(1+alpha,psi.*(beta+sl.*z1))./igamma(alpha,0));
ss2=exp(-N*sl.*z1).*(beta+sl.*z1).^(-1-alpha).*ss.* (beta.^alpha).*sl;


ss3=ss2+ sl.*(N+psi).*exp(-(N+psi).*sl.*z1) .*igamma(alpha,psi.*beta)/igamma(alpha,0) ;
y=exp(z).* ss3;


end