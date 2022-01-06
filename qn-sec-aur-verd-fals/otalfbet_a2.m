% --------------------------------------------------------------%
% Cálculo dos parâmetros alfa e beta  
  function [limalfa,limbeta,alfa,beta] = otalfbet_a2(xk,dxk,k,Xmax,Xmin,MAXITER)
% limalfa: Limite inferior para o valor de alfa
% limbeta: Limite superior para o valor de beta
% alfa   : Valor de alfa
% beta   : Valor de beta
% xk     : current soluction
% dxk    : direcao de busca
% Xmin   : Limite inferior para o xk
% Xmax   : Limite superior para o xk
% MAXITER: Numero maximo de iteracoes
% fatk   : parameter [0,1]
% --------------------------------------
% Definição do valor nulo.
  epslon = 1e-12;  
  fatk = 1;
  dxmax = Xmax-Xmin;
  s = fatk*(0.7 - 0.5*k/MAXITER);
  alfa = xk - s*dxmax;
  beta = xk + s*dxmax;
  sk = 0.7 + 0.2*k/MAXITER;
  for ii=1:length(Xmin)
      if alfa(ii,1) <= Xmin(ii,1)
         alfa(ii,1) = Xmin(ii,1) + (xk(ii,1)- Xmin(ii,1))*sk;
      end
      if beta(ii,1) >= Xmax(ii,1)
         beta(ii,1) = Xmax(ii,1) - (Xmax(ii,1)- xk(ii,1))*sk;
      end
  end

% ------------------------------------
% Gerar limite superior e inferior para alpha
  for ii=1:length(Xmax)
     if abs(dxk(ii,1)) > epslon
        k1(ii,1) = (alfa(ii,1) - xk(ii,1))/dxk(ii,1);
        k2(ii,1) = (beta(ii,1) - xk(ii,1))/dxk(ii,1);
     else
        k1(ii,1) = (alfa(ii,1) - xk(ii,1))/(sign(dxk(ii,1))*1.0e-6);
        k2(ii,1) = (beta(ii,1) - xk(ii,1))/(sign(dxk(ii,1))*1.0e-6);
     end
     if k1(ii,1) > k2(ii,1)
        aux1=k2(ii,1);
        k2(ii,1)=k1(ii,1);
        k1(ii,1)=aux1;
     end 
  end
  limalfa = max(k1);
  limbeta = min(k2);
  if limalfa > limbeta
     aux1=limbeta;
     limbeta=limalfa;
     limalfa=aux1;
  end
% Evitar que o valor de limalfa seja negativo.
  if limalfa < 0
     limalfa = 0;
  end 
% --------------------------------------------
% Deletar as variaveis que não são mais úteis.
  clear k1;
  clear k2;
  clear ftaux1;
  clear aux1;
  clear fdelta;
% ----------------------------------------