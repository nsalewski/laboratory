function [x,err,v]=myNewton(f,df,x0)
iteration=0;
tempx(1)=x0;
err=1;%Fehler wird auf einen beliebigen Wert größer als 10^(-12) gesetzt, damit die schleife startet
while (abs(err)>=(10^(-12))&&  iteration<50)
    iteration=iteration+1;
    k=iteration+1;
    tempx(k)=tempx(k-1)-(f(tempx(k-1))/(df(tempx(k-1))));
    err=abs(tempx(k)-tempx(k-1));
end
%Belegung der Rückgabevektoren erfolgt erst hier, da es so einfacher ist,
%einen Vektor angepasster Länge zu erzeugen.

x(1:iteration)=tempx(2:iteration+1);
err(1:iteration)=abs(tempx(2:iteration+1)-tempx(1:iteration));
v=f(x);
endfunction
