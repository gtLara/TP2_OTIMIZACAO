% ------------------------------------------------------------------------
% Calcula uma aproximação para a hessiana ou a inversa da hessiana (MMV).
  function IHk = otihmmv_X(fold,fnew,gkold,gknew,xkold,xknew,IHk,IIHES,Imetqn)
% ----------------------------------
% iteration - iteração corrente
% xkold     - vetor para cálculo de gkold.
% xknew     - vetor para cálculo de gknew.
% IHk       - entrada última estimativa da inversa da hessiana e
%             saída nova estimativa da inversa da hessiana.
% gkold     - gradiente em xkold.
% gknew     - gradiente em xknew.
% IIHES  : Define a operação a ser realizada.
%      1 : Calcula uma aproximação para a hessiana.
%     -1 : Calcula uma aproximação para a inversa da hessiana.
% -----------------------------------------
% Imetqn : Controla a escolha do algoritmo DFP ou BFGS.
%      0 : Algoritmo DFP (tetha = 0.0 , pho = 1.0).
%      1 : Algoritmo BFGS(tetha = 1.0 , pho = 1.0).
%      2 : Huang.
%      3 : Biggs.
% -----------------------------------------
% Interface de saída.
% IHk    : saída nova estimativa da hessiana ou da inversa da hessiana.
% -------------------------
% Precisao
  epslon = 1e-15;
% -----------------------------------------------------
% Icorr : Metodo Quase-Newton para funcoes quadraticas.
%     0 : DFP
%     1 : BFGS
%     2 : Biggs - Corrigir funcoes nao quadraticas
%     4 : Huang - Corrigir funcoes nao quadraticas.
% ------------------------------------------------
  Icorr = 0;
  tetha = 1.0; % Se BFGS
  pho   = 1.0;
  if Imetqn == 0 % Se DFP
     tetha = 0.0;
     pho = 1.0;
  elseif Imetqn == 2
     Icorr = 1;
  elseif Imetqn==3
     Icorr = 2;
  end
% ------------------------------------
  
% Verificar problema de divisão por zero.
  vetp = xknew - xkold;
  vety = gknew - gkold;

% ==========================================================
% Cálculo do fator eta para corrigir funções não qradráticas. 
% Huang
  if Icorr == 1
     tet=6*(fold-fnew)+3*(gkold+gknew)'*vetp;    
     vety=(1+tet/((vetp')*vety))*vety;
  end

% Biggs
  eta = 1;
  if Icorr == 2
     eps = (fnew-fold)/((vetp')*gkold);
     beta = (vetp')*gknew/((vetp')*gkold);
     if abs(2*beta + 1 -3*eps) > epslon
        eta = 0.5*(beta - 1)/(2*beta + 1 - 3*eps);
     end
  end
  
% Calcula uma aproximação para a matriz inversa da Hessiana.
  if IIHES == -1
     sig = (vetp')*vety;
     tau = (vety')*IHk*vety;
     if abs(sig) <= epslon || abs(tau) <= epslon
      % IHk =  eye(length(xknew));
        return
     end
     IHk = (IHk - IHk*vety*(vety')*IHk/tau +   tetha*(tau*vetp*(vetp')/sig^2 - ...
            IHk*vety*(vetp')/sig - vetp*((IHk*vety)')/sig + ...
            IHk*vety*(vety')*IHk/tau))*pho + ...
            eta*vetp*(vetp')/sig;        
  end      

% Calcula uma aproximação para a matriz Hessiana.
  if IIHES == 1
     sig = (vety')*vetp;
     tau = (vetp')*IHk*vetp;
     if abs(sig) <= epslon || abs(tau) <= epslon
      % IHk =  eye(length(xknew));
        return
     end
     IHk = (IHk - IHk*vetp*(vetp')*IHk/tau + tetha*(tau*vety*(vety')/sig^2 - ...
            IHk*vetp*(vety')/sig - vety*(IHk*vetp)'/sig + ...
            IHk*vetp*(vetp')*IHk/tau))*pho + ...
            eta*vety*(vety')/sig;
  end
% ----------------------------------------------------------------
% Verifica se matriz é positivo definida ou não (Decomp.CHOLESKY).
  [IER]= otdchol(IHk,epslon);
  if IER ~= 0
     disp(' Problema   Decomp. Cholesky  ');
   % Calcula uma aproximação para a matriz Hessiana ou inversa da Hessiana.
     IHk = eye(length(xknew));
  end
% =====================================================================
% =====================================================================
   function [IER]= otdchol(MATRIC,EPSLON)
% ---------------------------------------------------------
%  Verification possibilite Decomposition de CHOLESKY,
%  determination matrice d'entree definie positive ou non.
%  Interface d'entree
%  MATRIC : Matrice d'entree.
%  EPSLON : A specified tolerance (abs(X) < EPSLON). 
%  ---------------------------------------------
%  Interface de sortie.
%  IER   : Verification possibilite Decomposition de CHOLESKY.
%          0 : Decomposition de CHOLESKY possible,
%              Matrice d'entree definie positive.
%          1 : Decomposition de CHOLESKY impossible,
%              Matrice d'entree non definie positive.
%  ---------------------------------------
%  Possibilite Decomposition de CHOLESKY.
%  Matrice d'entree definie positive ou non.
   IER = 0;
%  --------------------------------------------------
%  Definição do número de linhas e colunas da matriz.
   [LINHA,COLUNA]=size(MATRIC);
   if (LINHA ~= COLUNA)
      disp(' ---------------------------------------');
      disp('  Erro fatal : Matriz não quadrada. ');
      disp('  Tecle <ctrl>C para terminar a execução. ');
      disp(' ---------------------------------------');
      disp('    '),
      pause;
      return;
   end
%  --------------------------------------
%  Verifica limite para parametro EPSLON.
   EPSLON1 = EPSLON;
%  ------------------------------------------------
%  Altera valores da diagonal da matriz de entrada.
   for I = 1 : LINHA
       MATRIC(I,I) = (1.0 - EPSLON1)*MATRIC(I,I);
   end
   for K = 1:LINHA 
       if (MATRIC(K,K) <= EPSLON1 )
        % Matriz de entrada não definida positiva.
          IER = 1;
          return;
       end
       MATRIC(K,K) = sqrt(MATRIC(K,K));
       for I = K+1:LINHA
           MATRIC(I,K) = MATRIC(I,K)/MATRIC(K,K);
           for J = K+1:I
               MATRIC(I,J) = MATRIC(I,J) - MATRIC(I,K)*MATRIC(J,K);
           end
       end
    end
% =====================================================================
% =====================================================================