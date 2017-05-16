x_ls = linspace(-1,1);%erzeugt die x-Werte, bezüglich derer die Funktion und die Approximationspolynome für das Plotten ausgewertet werden. 
delta=zeros(1,101);% dient zur Berechnung des maximalen Werts max|f(xi)-p(xi)|; xi ist dabei Element von Delta
for i=0:100 %erzeugt Gitter für die Betrachtung des maximalen Fehlers zwischen Funktion und Approximationspolynom
  delta(i+1)=-1+2*i/100;
end

%Hilfsfunktionen 

%berechnet die y-Werte zu gegebenen Stützstellen, x-Werten und c_i Vorfaktoren
function [y] = local_newton_polynom(x_plot,x_supps,c) 
	m=length(x_supps);
  y = c(length(x_supps)); %y wird auf das letzte c_i gesetzt
	for k=m-1:-1:1
		y = c(k)+(x_plot-x_supps(k)).*y;%es wird wieder quasi rückwärts gerechnet. So enthält jeder Summand im y außer dem ersten Summanden den Faktor (x-x_0), jeder außer den ersten beiden den Faktor (x-x_1) etc...
    %rückwarts gerechnet erhält man somit die gesamte Summe für y, indem an das y des vorherigen Schritts der nächstkleinere Faktor der Form (x-x_i) multipliziert wird, und das nächstkleinere c_i addiert wird.
   end
endfunction

%Runge-Funktion
function [f]=local_runge_function(x)
  f=(1)./(1+25.*x.^2);
endfunction

%Plottet die entsprechenden Graphen und berechnet den maximalen Fehler zur Runge-Funktion
function []=local_plot_fig_and_max_error(n,x_ls, filename, figtitle,delta,legend_position) 
x=zeros(1,n+1);
t=zeros(1,n+1);
for i=1:n+1
  x(i)=-1+((2*(i-1))/n);  
  t(i)=cos((2*(i-1)+1)*pi/(2*n+2));
end
figure
p=plot(x_ls,local_newton_polynom(x_ls,x,myNewtonInterpol(x,local_runge_function(x))),"r-",x_ls,local_newton_polynom(x_ls,t,myNewtonInterpol(t,local_runge_function(t))),"b-",x_ls,local_runge_function(x_ls),"g-", x,y=0,"r*",t,y=0,"b*");
title(figtitle)
legend(p([1,2,3,4,19]),'p(x) mit aequidistanten Stuetzstellen','p(x) mit Tschebyscheff-Knoten','f(x)','aequidistante Stuetzstellen','Tschebyscheff-Knoten','Location',legend_position)%erzeugt Legende
xlabel('x')
ylabel('f(x)')
print(filename)
%Ausgabe des Maximums
disp(figtitle)
max_aequidistant=max(local_runge_function(delta)-local_newton_polynom(delta,x,myNewtonInterpol(x,local_runge_function(x))))
max_tschebyscheff=max(local_runge_function(delta)-local_newton_polynom(delta,t,myNewtonInterpol(t,local_runge_function(t))))
endfunction

%Berechnung für n=7
local_plot_fig_and_max_error(7,x_ls, 'PA2-1-N7.pdf', 'n=7',delta, 'northeast')

%Berechnung für n=12
local_plot_fig_and_max_error(12,x_ls, 'PA2-1-N12.pdf', 'n=12',delta, 'south')

%Berechnung für n=17
local_plot_fig_and_max_error(17,x_ls, 'PA2-1-N17.pdf', 'n=17',delta, 'south')

%{
Beobachtung:
das zu den Tschebyscheff-Knoten berechnete Polynom konvergiert um Null langsamer gegen die Runge-Funktion als das berechnete Newton-Polynom zu äquidistanten Stützstellen.
Allerdings zeigt sich, besonders für x-Werte abseits der Null, dass das Newton-Polynom zu den Tschebyscheff-Knoten die Runge-Funktion deutlich besser approximiert.
Für n=12, bzw. n=17 wird zudem ersichtlich, dass das Newton-Polynom zu äquidistanten Stützstellen für Werte außerhalb eines kleinen Intervalls um Null, nicht gegen die Runge-Funktion konvergiert, sondern stattdessen oszilliert.
Bei der Berechnung des maximalen Fehlers wird zudem ersichtlich, dass die Tschebyscheff-Knoten nur einen sehr kleinen Fehler in der Approximation erzeugen, der Fehler für äquidistante Stützstellen allerdings wird für größere n auch immer größer.



Terminalausgabe:
>> myNewtonInterpolTest
n=7
max_aequidistant =  0.247358606559315
max_tschebyscheff =  0.391740284590228
n=12
max_aequidistant =  3.60527445056120
max_tschebyscheff =  0.0650551678196997
n=17
max_aequidistant =  4.06355074332311
max_tschebyscheff =  0.0559074779118303
%}