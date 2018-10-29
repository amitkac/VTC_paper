function y_theo =NoiseIntCdf(Noise_var,num,lambda_11,psi,x)


y_theo=gamcdf(x-Noise_var,num,1/lambda_11)+...
    (1-gamcdf(x-Noise_var,num,1/lambda_11)).*heaviside(x-Noise_var-psi);

end