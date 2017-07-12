g=@(x)(1/6*x^3+1/3*x^2+1/6);
x0=0;
k=8;%Es wurde k zu 7.34, also 8 Schritten in der Theorie bestimmt.
function [x,err]=myFixpunktIter(x0,g,k)
  str=sprintf("Gewählter Startwert %d",x0);
  disp(str)
  tempx(1)=x0;
  for c=1:k
    tempx(c+1)=g(tempx(c));
  end
  err=abs(tempx(2:k+1)-tempx(1:k));
  x=tempx(2:k+1);
endfunction
[x,err]=myFixpunktIter(x0,g,k);
for i=1:k
  str=sprintf("x_%d=%.8f und |x_%d-x_%d|=%.12f",i,x(i),i,i-1,err(i));
  disp(str)
end
%{
Terminalausgabe:
Gewählter Startwert 0
x_1=0.16666667 und |x_1-x_0|=0.166666666667
x_2=0.17669753 und |x_2-x_1|=0.010030864198
x_3=0.17799348 und |x_3-x_2|=0.001295950505
x_4=0.17816708 und |x_4-x_3|=0.000173600503
x_5=0.17819044 und |x_5-x_4|=0.000023362545
x_6=0.17819359 und |x_6-x_5|=0.000003145993
x_7=0.17819401 und |x_7-x_6|=0.000000423674
x_8=0.17819407 und |x_8-x_7|=0.000000057057
%}
