% ======================================================
% ======================================================  
  function [lag,dlag] = fun_rosensuzuki_irr(xk)
% ------------------------------------------------------ 
% Saídas:
% fobj  : valor da função objetivo para a iteração.
% dfobj : gradiente da função objetivo.
% gj    : Valores das funções de restrição.
% mdh   : Matriz Jacobiabno para as restrições.
% gsinal: vetor define sinal >= 0 ou <= 0 para as restrições.
% ------------------------------------------------------
  
% RosenSuzukiProblem: Problema de minimização convexa restrita 
% Conforme J. B. Rosen e S. Suzuki; Comunications of the ACM, vol. 8
% n. 2, Fevereiro de 1965
% --------------------------------------------
% Entradas:
% x0    - solução ótima para o problema primal
% u0    - solução ótima para o problema dual
% --------------------------------------------
% Problema de Rosen-Suzuki
% Passo 1: escolha x0 dentro do domínio e ui >= 0
% parâmetros de entrada
  nxk = length(xk);
  
% define as soluções x0 (primal) e u0(dual) para gerar o problema de 
% Rosen-Suzuki
  for i1 =1:nxk-1
      x0(i1,1) = (-1)^i1;
      u0(i1,1) = (-1)^i1 + 1;
  end
  x0(nxk,1) = (-1)^nxk;
  
% -----------------------------------------------------------------------
% Passo 2: escolha bi tal que hi(x0) = 0 para ui > 0 e que hi(x0) < 0
% para uj = 0
% qj = x'*Qj*x + ai'*x + bj
  mdg = [];
  for j1 =1:nxk -1
      for i1 =1:nxk
          vx(i1,1) = 2 - (-1)^(i1+j1);
          ai(i1,1) = 1 + (-1)^j1 + (-1)^i1;
      end
      Qj = diag(vx);
      gj(j1,1) = x0'*Qj*x0 + ai'*x0;
      mdg = [mdg (2*Qj*x0 + ai)];
      if u0(j1,1) > 0
         bj(j1,1) = -gj(j1,1);     % calcula bi para gj(i) = 0
      else
         bj(j1,1) = -gj(j1,1) - 1; % soma -1 para fazer gj(i) < 0
      end
  end

% -----------------------------------------------------------------------
% Passo 3: calculo de c = -grad(theta(x0)) - somatório(ui*grad(qi(x0)))
% fobj =  x'*Q0*x + c'*x0
  for i1 =1:nxk
      vx(i1,1) = 2 - (-1)^i1;
  end
  Q0 = diag(vx);
  dfobj = 2*Q0*x0;
  c = - dfobj - sum(mdg*u0,2);
% ------------------------------------------------------
  mdg = [];
  for j1 =1:nxk -1
      for i1 =1:nxk
          vx(i1,1) = 2 - (-1)^(i1+j1);
          ai(i1,1) = 1 + (-1)^j1 + (-1)^i1;
      end
      Qj = diag(vx);
      gj(j1,1) = xk'*Qj*xk + ai'*xk + bj(j1,1);
      mdg = [mdg (2*Qj*xk + ai)];
      gsinal(j1,1) = -1;
  end 
  dfobj= 2*Q0*xk + c;
  fobj = xk'*Q0*xk + c'*xk;
% Lagrangean function
  lag = fobj + gj'*u0;
  dlag = dfobj + mdg*u0;
  end
% ======================================================
% ======================================================