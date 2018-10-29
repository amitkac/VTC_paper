function y_theo =NoiseIntPdf(Noise_var,num,lambda_11,psi,x)



y_theo=gampdf(x-Noise_var,num,1/lambda_11)+dirac(x-Noise_var-psi)- (heaviside(x-Noise_var-psi).*...
    gampdf(x-Noise_var,num,1/lambda_11))-gamcdf(x-Noise_var,num,1/lambda_11).*...
    dirac(x-Noise_var-psi);


end