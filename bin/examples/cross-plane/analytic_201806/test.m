clc;
Kn=1;
fun=@(mio) (mio.^2);% ( exp (mio)) ; %*(-1.^n) )./(1+(Kn(i)*m).^2*(n*pi).^2);
Integ_matlab=integral(fun,0,1)
Ngauss=4;
[PP WW]=GaussIntegrate(-1,1,Ngauss); %GaussIntegrationPoints
integ_num=0;
for i=1:Ngauss
    integ_num=integ_num+1/2*(PP(i)^2)*WW(i) ;%exp(PP(i))*WW(i) ; 
end
integ_num