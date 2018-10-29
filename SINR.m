function y_theo=SINR(Noise_var,num,lambda_11,sl,psi,z)

alpha=num;
  beta=lambda_11;
  % z=0:0.1:10;
 % sl=lambda_pu/p_pu;
  N=Noise_var;
%  s1=igamma(1+alpha,(alpha-N)  
%k=beta+sl.*z; 
q=psi;


ss= alpha+ N.*(beta+sl.*z).*(1-igamma(alpha,psi.*(beta+sl.*z))...
       ./igamma(alpha,0))-(igamma(1+alpha,psi.*(beta+sl.*z))./igamma(alpha,0));
ss2=exp(-N*sl.*z).*(beta+sl.*z).^(-1-alpha).*ss.* (beta.^alpha).*sl;


y_theo=ss2+ sl.*(N+psi).*exp(-(N+psi).*sl.*z) .*igamma(alpha,q.*beta)/igamma(alpha,0) ;

end