function [x,err,v]=mybisect(f,a,b)
  %die Variable e wurde auf err geändert, da e in Octave für die Exponetialfunktion reserviert ist.
assert(f(a)*f(b)<0, 'Nach dem Zwischenwertsatz befindet sich im Intervall keine Nullstelle. Berechnung abgebrochen')
%assert wirft error, wenn condition falsch.
iteration=0;
err=1;%Beliebiger startwert größer 10^(-12)
temx(1)=(a+b)/2;
while(iteration<9999&&abs(err)>=10^(-12))
  iteration=iteration+1;
  if(f(temx(iteration))==0)
    break;
  elseif(f(a)*f(temx(iteration))>0)
    %Wenn f(a) und f(intervallmitte) das gleiche Vorzeichen haben, ihr Produkt also positiv ist,
    %dann kann die Hälfte des Intervalls "abgeschnitten" werden, da sich die Nullstelle hier nicht befinden wird.
    a=temx(iteration);%a auf alte Intervallmitte setzen.
    temx(iteration+1)=(a+b)/2;%neue Intervallmitte bestimmen
  else
    b=temx(iteration);%Ist das Vorzeichen in elseif negativ, kann die andere Hälfte des Intervalls abgeschnitten werden.
    temx(iteration+1)=(a+b)/2;
  end
  err=abs(temx(iteration+1)-temx(iteration));%Berechnung des Fehlers. In iteration ist hierbei k-1 gespeichert.
end
%Rückgabe aller berechneter Werte
x(1:iteration)=temx(2:iteration+1);
err(1:iteration)=abs(temx(2:iteration+1)-temx(1:iteration));
v=f(x);
endfunction
