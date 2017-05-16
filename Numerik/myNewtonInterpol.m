## Author: julian <julian@julian-ThinkPad-L450>
## Created: 2017-05-12
%MYNEWTONINTERPOL Berechnet die Koeffizienten c_i des Newton-Interpolationspolynom
% MYNEWTONINTERPOL Für gegebene Stützstellen x_i und den zugehörigen Funktionswerten f_i werden die Koeffizienten c_i des Newton-Interpolationspolynom berechnet. Verwendung: c_i=myNewtonInterpol(x_i,f_i).
function [c] = myNewtonInterpol(x,f)
  c = f; %Zunächst werden alle Funktionswerte in den c_i gespeichert
  n = length(x); %n wird auf die Anzahl der Stützstellen gesetzt
	for(l=2:n) %l läuft über die ganze Länge des Vektors, außer c(1), da c(1) ja nicht mehr geändert wird, und direkt c(1)=f(1) gilt.
		for(k=n:-1:l)% siehe Anmerkung unten
			c(k) = (c(k)-c(k-1)) / (x(k)-x(k-l+1));
		end
	end
endfunction
%Der c-Vektor wird nun quasi von unten nach oben veraendert. Der k-te Wert ergibt sich dabei jeweils aus der Differenz des urspruenglichen k-ten Wertes mit dem k-1-ten Wert
%dividiert durch die Differenz aus zugeordneten, k-ten x-Wert und dem k-l+1-ten x-Wert. Ebenso wird für den k-1-ten Wert verfahren, etc. In jedem Folgeschritt wird schließlich über das "l" ein Element des Vektors weniger veraendert.
%So ergibt sich sukzessize der c_i Vektor 
