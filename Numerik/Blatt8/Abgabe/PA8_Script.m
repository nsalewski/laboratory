format long
f=@(x)((cos(2*x)).^2-x.^2);
df=@(x)(-2*(sin(4*x))-2*x);
%Intervallgrenzen für Bisektion
a=0;
b=0.75;
[xb,err_bisection,vb]=mybisect(f,a,b);

[x,err_newton,v]=myNewton(f,df,b);%Es wird b für x0 verwendet, da diese den gleichen Wert haben.
err_newton
%Fehler werden halblogarithmisch geplottet
n_newton=(1:length(err_newton));
n_bisection=(1:length(err_bisection));
semilogy(n_newton,err_newton,n_bisection,err_bisection);
legend('Fehler Newton-Verfahren','Fehler Bisektionsverfahren')
xlabel('Anzahl n der Iterationsschritte');
ylabel('absoluter Fehler');
title('Fehler des Newton-Verfahrens und des Bisektionsverfahrens zur Bestimmung der Nullstelle \n in Abhaengigkeit der Iterationsschritte')
grid on
grid minor off
print('PA8.3.fig')
