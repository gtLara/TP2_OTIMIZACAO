% ==========================================
% Busca o valor do alfa otimo para a passada 
% utilizando o metodo da Se��o Aurea.
  function [alfaotim,icfunc] = otgoldsc_a22(funcao,icfunc,isa_FV,Hess,...
            XK,dxk,k,fold,fnew,dfold,dfnew,Xmin,Xmax,MAXITER,a)
% ---------------------------
% Par�metros de sa�da         
% alfaotim: Alfa �timo
% fim     : Flag para a rotina.
% XK      % Historico dos xk's
% dxk     : Dire��o de busca
% k       : Numero da itera��o
% fold    : Valor hist�rico da fun��o.
% fnew    : Valor atual da fun��o.
% dfold   : Valor hist�rico do gradiente da fun��o.
% dfnew   : Valor atual do gradiente da fun��o.
% Xmin    : Limite inferior para o intervalo do alfa
% Xmax    : Limite superior para o intervalo do alfa
% MAXITER : Numero maximo de iteracoes para o algoritmo.
% isa_FV = 0 => Falsa se��o aurea.
% isa_FV = 1 => Verdadeira se��o aurea.

% Inicializa variaveis
  epslon=1e-14; % precisao para a secao aurea

% C�lculo dos par�metros alfa e beta
  [limalfa,limbeta,alfa,beta] = otalfbet_a2(XK(:,size(XK,2)),dxk,k,Xmax,Xmin,MAXITER);
  limalfa = 0;
  alfasup = limalfa + 0.618*(limbeta-limalfa);
  alfainf = limbeta - 0.618*(limbeta-limalfa);

% -----------------------------
% Estoca o valor inicial de xk.
  xk = XK(:,size(XK,2));
  xk1= XK(:,size(XK,2)-1);
  xk2= xk;

% Flag para definir numero finito de loops para a Se��o Aurea.
  nloops = 100;
  fim = 0;
  
% ---------------------------
% Gera valor de xk em alfasup.
  xksup = xk2 + alfasup*dxk;
% Se isa_FV igual zero Falsa se��o aurea
% Se isa_FV igual um verdadeira se��o aurea;
  if isa_FV == 0
   % Funcao gera aproximacoes para as funcoes para o metodo da Se��o Aurea.
     [faprx2] = otapp_a2(Hess,XK,k,xksup,fnew,dfnew);
  else
    % ==========================================
    % ==========================================
    % -------------------------------------------------
    % fob   : valor da fun��o objetivo para a itera��o.
    % dfobj : gradiente da fun��o objetivo.
    % fun��o: fun��o a ser avaliada
    % xk    : Ponto onde se deseja avaliar a fun��o
      [faprx2,dfobjnew] = feval(funcao,xksup,a);
      icfunc = icfunc + 1;
    % ==========================================
    % ==========================================
  end
% ---------------------------    
% Gera valor de xk em alfainf.
  xkinf = xk2 + alfainf*dxk;
   if isa_FV ==0
    % Funcao gera aproximacoes para as funcoes para o metodo da Se��o Aurea.
      [faprx1] = otapp_a2(Hess,XK,k,xkinf,fnew,dfnew);
   else
      % ==========================================
      % ==========================================
      % -------------------------------------------------
      % fob   : valor da fun��o objetivo para a itera��o.
      % dfobj : gradiente da fun��o objetivo.
      % fun��o: fun��o a ser avaliada
      % xk    : Ponto onde se deseja avaliar a fun��o
      [faprx1,dfobjnew] = feval(funcao,xkinf,a);
      icfunc = icfunc + 1;
      % ==========================================
      % ==========================================
  end
  while (0.618*(limbeta-limalfa)) > epslon & fim <= nloops
        if faprx1 > faprx2 
           limalfa = alfainf;
           alfainf = alfasup;
           alfasup = limalfa + 0.618*(limbeta-limalfa);
           faprx1 = faprx2;
        
         % ----------------------------
         % Gera valor de xk em alfasup.
           xksup = xk2 + alfasup*dxk;                     
          
          if isa_FV == 0
              % Funcao gera aproximacoes para as funcoes para o metodo da Se��o Aurea.
              [faprx2] = otapp_a2(Hess,XK,k,xksup,fnew,dfnew);
          else
              % ==========================================
              % ==========================================
              % -------------------------------------------------
              % fob   : valor da fun��o objetivo para a itera��o.
              % dfobj : gradiente da fun��o objetivo.
              % fun��o: fun��o a ser avaliada
              % xk    : Ponto onde se deseja avaliar a fun��o
              [faprx2,dfobjnew] = feval(funcao,xksup,a);
               icfunc = icfunc + 1;
              % ==========================================
              % ==========================================
          end           
           
         % Flag para impedir numero infinito de loops
           fim = fim + 1;
         else 
            limbeta = alfasup;
            alfasup = alfainf;
            alfainf = limbeta - 0.618*(limbeta-limalfa);
            faprx2 = faprx1;
            
          % ---------------------------
          % Gera valor de xk em alfainf.
            xkinf = xk2 + alfainf*dxk;
            
            if isa_FV ==0
                % Funcao gera aproximacoes para as funcoes para o metodo da Se��o Aurea.
                [faprx1] = otapp_a2(Hess,XK,k,xkinf,fnew,dfnew);
            else
                % ==========================================
                % ==========================================
                % -------------------------------------------------
                % fob   : valor da fun��o objetivo para a itera��o.
                % dfobj : gradiente da fun��o objetivo.
                % fun��o: fun��o a ser avaliada
                % xk    : Ponto onde se deseja avaliar a fun��o
                [faprx1,dfobjnew] = feval(funcao,xkinf,a);
                 icfunc = icfunc + 1;
                % ==========================================
                % ==========================================
            end
           
          % Flag para impedir numero infinito de loops
            fim = fim + 1;
         end
  end
% ------------------------------------------                  
% Valor final para o alfa �timo
  if faprx1 > faprx2
     limalfa = alfainf;
  else
     limbeta = alfasup;
  end
  alfaotim = abs(limalfa+limbeta)/2;
% ================================================================
% Gera aproximacao quadratica para a funcao.
  function [funaprox] = otapp_a2(Hess,XK,k,xkdavez,fnew,dfnew)
% Par�metros de sa�da         
% fim     : Flag para a rotina.
% XK      % Historico dos xk's
% k       : Numero da itera��o
% fnew    : Valor atual da fun��o.
% dfnew   : Valor atual do gradiente da fun��o.
% xkdavez : Valor atual de xk.
% -------------------------------------
% Gera uma aproxima��o para a fun�ao em xkdavez.
  dx = xkdavez - XK(:,k);
  funaprox = fnew + dfnew'*dx + 0.5*dx'*Hess*dx;
% =====================================================================