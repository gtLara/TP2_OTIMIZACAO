 function [fobj,dfobj]=fex3(x, a)
% funcao a ser minimizada (2 variáveis)

% ===========================================================
% Função objetivo.
  fobj = -8*x(1,:).*x(2,:) + (4/a^2)*x(2,:).*x(1,:).^3  + (4/a^2)*x(1,:).*x(2,:).^3;
  
% Gradiente da função objetivo.
  dfobj(1,1) = -8*x(2,1) + (12/a^2)*x(2,1)*x(1,1)^2 + (4/a^2)*x(2,1)^3;
  dfobj(2,1) = -8*x(1,1) + (12/a^2)*x(1,1)*x(2,1)^2 + (4/a^2)*x(1,1)^3;

% ----------------------
