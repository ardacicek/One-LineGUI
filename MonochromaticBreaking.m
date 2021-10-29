function [H0,L0,T,db, Hb,Lb, alfab]=MonochromaticBreaking(H0,T,alfa0,m)
alfa0abs = abs(deg2rad(alfa0));
L0=1.56*T^2;
dbassumed=H0;
Lassumed=L0;
db=0;
a=43.8*(1-exp(-19*m));
b=1.56/(1+exp(-19.5*m));
errorL=1;
errordb=1;
c=0;
while errordb>0.000000001
    while errorL>0.0001
        Lb=L0*tanh(2*pi*dbassumed/Lassumed);
        errorL=abs(Lassumed-Lb);
        Lassumed=Lb;
    end
   alfab=asin(tanh(2*pi*dbassumed/Lb)*sin(alfa0abs));
   Krb=sqrt(cos(alfa0abs)/cos(alfab));
   H0prime=H0*Krb;
   Hb=0.53*H0prime*(H0prime/L0)^(-0.24);
   db=Hb/(b-a*Hb/(9.81*T^2));
   errordb=abs(dbassumed-db);
   dbassumed=db;
end
if alfa0 < 0
    alfab = -alfab;
end