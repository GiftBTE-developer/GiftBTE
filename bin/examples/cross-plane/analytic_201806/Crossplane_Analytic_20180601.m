clc
%In the name of Allah, the compassionate the most merciful;
%Analytic solution of cross-plane problem-1D
%f_N(x_hat)=f0/2+sum(fm *cos(m*PI*x_hat)
format short ;
%double(x)
L=5E-9  %1.06812E-5
%Ngauss=4; 
%L=1.068E-8 ,-6 ,-5: k
%=8.52 % 115.0314 , 151.8587 //  5bandsiso
%[PP WW]=GaussIntegrate(-1,1,Ngauss); %GaussIntegrationPoints
%L=1.068E-7 ,-6 ,-5: k=43.6357 , 104.9934,147.2315//  50bands
filename = 'Input-dispersion-relation-fp-15.dat'; %5-isotropic.dat' ; %-20-ProfYi.dat';
[dispertioncurve,delimiterOut]=importdata(filename);
vg  =dispertioncurve(:,1)  ;
Tau =dispertioncurve(:,2)  ;
C=dispertioncurve(:,3)     ;
Lr  =dispertioncurve(:,4)  ;
Kn  =vg.*Tau/L             ;
Nband=length(Kn)           
% Kn(1)=0.01;
% C(1)=1    ;
% Tau(1)=1  ;
% Kn(2)=100 ;
% C(2)=1    ;
% Tau(2)=1  ;  
DeltaT=1    ;
Mm=100      ;% Mm=Nn
Nn=100      ;% Nn=100: 138.7537
Nt=100; %number of domain division points
%integrate of C/Tau
sum_C_Tau=0;
for ib=1:Nband
    sum_C_Tau=sum_C_Tau+C(ib)/Tau(ib);
end

f=zeros(Mm,1); 
%f0
sum_CKn_Tau=0;
for ib=1:Nband
    fun=@(mio) mio.*( exp (-1./(Kn(ib)*mio))); 
    E3=integral(fun,0,1);
    sum_CKn_Tau=sum_CKn_Tau+C(ib)/Tau(ib)*Kn(ib)*DeltaT*(1-2*E3);
end 
f0=sum_CKn_Tau/(2*sum_C_Tau);
f(1)=f0;

%fm
for m=2:Mm
    %m=m-1;
    sum_CKnm_Tau=0;
    for ib=1:Nband
        fun=@(mio)  mio.*( 1-(exp (-1./(Kn(ib)*mio)))*((-1).^(m-1)))./ (1+(Kn(ib)*mio).^2*((m-1)*pi).^2);
        integ=integral(fun,0,1);
        sum_CKnm_Tau=sum_CKnm_Tau+C(ib)/Tau(ib)*Kn(ib)*DeltaT*integ; 
    end
    f(m)=sum_CKnm_Tau/sum_C_Tau; %fm
end

k=zeros(Nn,Nn); 
%k00
sum_CKn_Tau2=0;
for ib=1:Nband
    fun=@(mio)  mio.*( exp (-1./(Kn(ib)*mio))) ; %*(-1.^n) )./(1+(Kn(i)*m).^2*(n*pi).^2);
    E3=integral(fun,0,1);
    sum_CKn_Tau2=sum_CKn_Tau2+C(ib)/Tau(ib)*Kn(ib)*(2/Kn(ib)-1+2*E3); 
end
k(1,1)=2*sum_CKn_Tau2/sum_C_Tau; %k00

%kn0
for n=2:Nn
    %n=n-1;
    sum_CKnm_Tau2=0;
    for ib=1:Nband
        fun=@(mio)  mio.* (exp (-1./(Kn(ib)*mio))-1)./ (1+(((Kn(ib)*mio).^2).*(((n-1) *pi).^2)));
        integ=integral(fun,0,1);
        sum_CKnm_Tau2=sum_CKnm_Tau2+C(ib)/Tau(ib)*Kn(ib)*(((-1)^(n-1))+1)*integ; 
    end
    k(n,1)=2*sum_CKnm_Tau2/sum_C_Tau; %km0
end

%k0n
for n=2:Nn
    %n=n-1;
    sum_CKnm_Tau3=0;
    for ib=1:Nband
        fun=@(mio)  mio.* (exp (-1./(Kn(ib)*mio))-1)./ (1+((Kn(ib)*mio).^2) .*(((n-1) *pi).^2));
        integ=integral(fun,0,1);
        sum_CKnm_Tau3=sum_CKnm_Tau3+C(ib)/Tau(ib)*Kn(ib)*(((-1)^(n-1))+1)*integ ;
    end
    k(1,n)=2*sum_CKnm_Tau3/sum_C_Tau; %k0n
end

%kmn
for m=2:Mm
    for n=2:Nn
        if (m~=n)
            %n=n-1;
            %m=m-1;
            sum_CKnm_Tau4=0;
            for ib=1:Nband
                fun=@(mio)  mio.* (exp (-1./(Kn(ib)*mio))*(( (-1)^(n-1) + (-1)^(m-1) ))-(1+( (-1)^((n-1)+(m-1)) )))./(  (1+((Kn(ib)*mio).^2 )*(((n-1)*pi).^2)).*(1+((Kn(ib)*mio).^2 )*(((m-1) *pi).^2)) );
                integ=integral(fun,0,1);
                sum_CKnm_Tau4=sum_CKnm_Tau4+C(ib)/Tau(ib)*Kn(ib)*integ; 
            end
            k(m,n)=2*sum_CKnm_Tau4/sum_C_Tau;
        end
    end
end    
    
%knn n~=0 :  n not zero
for n=2:Nn 
    sum_C_Tau2=0;
    for ib=1:Nband
        fun=@(mio)  mio.*(( (exp(-1./(Kn(ib).*mio)) * ((-1).^(n-1))) -1) ./ (  (1+((Kn(ib)*mio).^2) .*(((n-1) *pi).^2) ).^2  )) ;  %mio.* (exp (-1./(Kn(ib)*mio))*(( (-1)^(n-1) + (-1)^(m-1) ))-(1+( (-1)^((n-1)+(m-1)) )))./(  (1+((Kn(ib)*mio).^2 )*(((n-1)*pi).^2)).*(1+((Kn(ib)*mio).^2 )*(((m-1) *pi).^2)) );
        integ=integral(fun,0,1);
        sum_C_Tau2=sum_C_Tau2+C(ib)/Tau(ib) *( (atan((n-1)*pi*Kn(ib))/((n-1)*pi*Kn(ib))) + (2*Kn(ib)*integ) ); 
    end
    k(n,n)=2*sum_C_Tau2/sum_C_Tau;
end

%F matrix
F=f;

%A matrix
A=zeros(Nn,Nn);
A(1,1)=1-k(1,1)/4;
for n=2:Nn
    A(1,n)=-k(1,n)/2;
end

for n=2:Nn
    A(n,1)=-k(n,1)/4;
end

for n=2:Nn
    for m=2:Nn   
        if (n==m)
           A(n,m)=1-k(n,n)/2;
        else
           A(n,m)=-k(n,m)/2;     
        end   
    end
end
%AX=F
A_1=(A^-1);
F;
X=(A_1)*F;


x_h(1)=0; %T(1) ?
x_h(Nt+1)=0; %T(11)?
L0=1;
for nt=1:Nt+1
    x_h(nt)=(nt-1)*L0/Nt;
end
T(1:Nt+1)=0;
disp('Temperature')
for nt=1:Nt+1
    for n=1:Nn 
        T(nt)=T(nt)+X(n)*cos((n-1)*pi*x_h(nt));
    end
    T(nt)=T(nt)-0.5;
    %Prnt = sprintf(' Temp   %d   %d', nt, T(nt));
    Pr = ['x: ', num2str(x_h(nt)),'      ', num2str(T(nt))];
    disp(Pr)
end

%Suppression calculation
SF=zeros(Nband);
for ib=1:Nband
    fun=@(mio)  (mio.^2).*( exp (-1./(Kn(ib)*mio))) ; 
    E4=integral(fun,0,1);     
    sum_xm_mio=0;
    for n=2:Nn
        fun=@(mio)  ((mio.^2).*(  1+(exp (-1./(Kn(ib)*mio)))  )./ (1+(((Kn(ib)*mio).^2)*(((n-1)*pi).^2))));
        integ=integral(fun,0,1);
        sum_xm_mio=sum_xm_mio+X(n)*(1-(-1)^(n-1))*integ;
    end
    SF(ib)=(1/2)-(3/2)*E4+(3/2)*sum_xm_mio;
   % Prsf = ['iband= ', num2str(ib),'  , k= ', num2str(SF(ib))]
   % disp(Prsf)
end
ks=0;
for ib=1:Nband
ks=ks+1/3*C(ib)*(vg(ib)^2)*Tau(ib)*SF(ib) ;  
hfband=1/3*C(ib)*(vg(ib)^2)*Tau(ib)*SF(ib);    
Prshf = ['iband= ', num2str(ib),'  , k= ',num2str(hfband) ];
    disp(Prshf)
end 
ks
data=[x_h',T'];
save('T1000nm','data');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %Yi Formula
% k2=0;
% for ib=1:Nband
%     sum_T_Kn=0;
%     fun=@(mio)  ( exp (-1./(Kn(ib)*mio))) ; 
%     E2=integral(fun,0,1);  
%     for ng=1:Ngauss     
%         for nt=1:Nt
%             if x_h(nt)<PP(ng)<x_h(nt+1)
%                ntemp=nt;
%             end    
%            % Pr = ['x: ', num2str(nt),'      ', num2str(Temp)];
%            % disp(Pr)
%         end
%         Temp=T(ntemp)-0.5
%         sum_T_Kn=sum_T_Kn+1/2*(integ/Kn(ib)*Temp)*WW(ng);
%     end
%     
%     k2=k2+L/(4*DeltaT)*(C(ib)*vg(ib))*(1-2*sum_T_Kn) ;
% end 
% k2




% for iband=1:Nband
%    fun = @(x,Y) ((1.-x.^2.).*(2.-exp(-Y./(Kn(iband)*x))-exp(-(1.-Y)./(Kn(iband)*x))));
%     q = integral2(fun,0,1,0,1);
% end
% integ=integ+3/4*Kr(iband)*q;     
% krc=integ;
% X = sprintf('keff/kbulk [%d ]=%d.',Kn(1),krc);
% disp (X);