  function [fobj,dfobj]=fexlivro(x, a)
% funcao a ser minimizada (2 variáveis)

% ===========================================================
% Limite superior do intervalo de busca 
  % x_sup = [10; 10]; 
  % a_sup = max(roots([36*x_sup(1,1)*x_sup(2,1) 144*x_sup(2,1)+48*x_sup(1,1) 48]));
     
% % Limite inferior do intervalo de busca
  % x_inf = [-10; -10]; 
  % a_inf = min(roots([36*x_inf(1,1)*x_inf(2,1) 144*x_inf(2,1)+48*x_inf(1,1) 48]));
  % Val_a = 0;
  % if Val_a == 0
  %  % Definição do valor a = 0
  %    a = (a_sup + a_inf)/2;
  % end

  % if Val_a == 1
  %  % Definição do valor a positivo próximo de zero.
  %    a = (a_sup + 0.95*a_inf)/2;
  % end
  % if Val_a == -1

  %  % Definição do valor a negativo próximo de zero.
  %    a = (0.95*a_sup + a_inf)/2;
  % end
%    xstar(1,1) = min(roots([0.1875*a^3 3*a^2 (0.25*a^2+14*a) (2*a+4) ((1/12)*a+8/6)]));
%    xstar(2,1) = min(roots([0.1875*a^3 a^2 ((4/3)*a+6*a) 4 2]));
% ===========================================================
% Função objetivo.
  fobj = 12*x(1,1)^2 + 4*x(2,1)^2 - 12*x(1,1)* x(2,1) + 2*x(1,1)+ a*(x(1,1)^3+x(2,1)^3);
  
% Gradiente da função objetivo.
  dfobj(1,1) = 24*x(1,1) - 12*x(2,1) + 2 + 3*a*x(1,1)^2;
  dfobj(2,1) = 8*x(2,1)  - 12*x(1,1) + 3*a*x(2,1)^2;
% ----------------------
