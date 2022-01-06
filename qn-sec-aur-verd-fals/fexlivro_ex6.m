function [fobj,dfobj]=fexlivro_ex6(x)
% funcao a ser minimizada (2 variáveis)

global a

%    xstar(1,1) = min(roots([0.1875*a^3 3*a^2 (0.25*a^2+14*a) (2*a+4) ((1/12)*a+8/6)]));
%    xstar(2,1) = min(roots([0.1875*a^3 a^2 ((4/3)*a+6*a) 4 2]));
% ===========================================================
% Função objetivo.
  fobj = 12*x(1,1)^2 + 4*x(2,1)^2 - 12*x(1,1)* x(2,1) + 2*x(1,1)+ a*(x(1,1)^3+x(2,1)^3);
  
% Gradiente da função objetivo.
  dfobj(1,1) = 24*x(1,1) - 12*x(2,1) + 2 + 3*a*x(1,1)^2;
  dfobj(2,1) = 8*x(2,1)  - 12*x(1,1) + 3*a*x(2,1)^2;
% ----------------------