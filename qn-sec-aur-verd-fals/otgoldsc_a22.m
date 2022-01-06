% ==========================================
% Busca o valor do alfa otimo para a passada 
% utilizando o metodo da Seção Aurea.
  function [alfaotim,icfunc] = otgoldsc_a22(funcao,icfunc,isa_FV,Hess,...
            XK,dxk,k,fold,fnew,dfold,dfnew,Xmin,Xmax,MAXITER)
% ---------------------------
% Parâmetros de saída         
% alfaotim: Alfa ótimo
% fim     : Flag para a rotina.
% XK      % Historico dos xk's
% dxk     : Direção de busca
% k       : Numero da iteração
% fold    : Valor histórico da função.
% fnew    : Valor atual da função.
% dfold   : Valor histórico do gradiente da função.
% dfnew   : Valor atual do gradiente da função.
% Xmin    : Limite inferior para o intervalo do alfa
% Xmax    : Limite superior para o intervalo do alfa
% MAXITER : Numero maximo de iteracoes para o algoritmo.
% isa_FV = 0 => Falsa seção aurea.
% isa_FV = 1 => Verdadeira seção aurea.

% Inicializa variaveis
  epslon=1e-14; % precisao para a secao aurea

% Cálculo dos parâmetros alfa e beta
  [limalfa,limbeta,alfa,beta] = otalfbet_a2(XK(:,size(XK,2)),dxk,k,Xmax,Xmin,MAXITER);
  limalfa = 0;
  alfasup = limalfa + 0.618*(limbeta-limalfa);
  alfainf = limbeta - 0.618*(limbeta-limalfa);

% -----------------------------
% Estoca o valor inicial de xk.
  xk = XK(:,size(XK,2));
  xk1= XK(:,size(XK,2)-1);
  xk2= xk;

% Flag para definir numero finito de loops para a Seção Aurea.
  nloops = 100;
  fim = 0;
  
% ---------------------------
% Gera valor de xk em alfasup.
  xksup = xk2 + alfasup*dxk;
% Se isa_FV igual zero Falsa seção aurea
% Se isa_FV igual um verdadeira seção aurea;
  if isa_FV == 0
   % Funcao gera aproximacoes para as funcoes para o metodo da Seção Aurea.
     [faprx2] = otapp_a2(Hess,XK,k,xksup,fnew,dfnew);
  else
    % ==========================================
    % ==========================================
    % -------------------------------------------------
    % fob   : valor da função objetivo para a iteração.
    % dfobj : gradiente da função objetivo.
    % função: função a ser avaliada
    % xk    : Ponto onde se deseja avaliar a função
      [faprx2,dfobjnew] = feval(funcao,xksup);
      icfunc = icfunc + 1;
    % ==========================================
    % ==========================================
  end
% ---------------------------    
% Gera valor de xk em alfainf.
  xkinf = xk2 + alfainf*dxk;
   if isa_FV ==0
    % Funcao gera aproximacoes para as funcoes para o metodo da Seção Aurea.
      [faprx1] = otapp_a2(Hess,XK,k,xkinf,fnew,dfnew);
   else
      % ==========================================
      % ==========================================
      % -------------------------------------------------
      % fob   : valor da função objetivo para a iteração.
      % dfobj : gradiente da função objetivo.
      % função: função a ser avaliada
      % xk    : Ponto onde se deseja avaliar a função
      [faprx1,dfobjnew] = feval(funcao,xkinf);
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
              % Funcao gera aproximacoes para as funcoes para o metodo da Seção Aurea.
              [faprx2] = otapp_a2(Hess,XK,k,xksup,fnew,dfnew);
          else
              % ==========================================
              % ==========================================
              % -------------------------------------------------
              % fob   : valor da função objetivo para a iteração.
              % dfobj : gradiente da função objetivo.
              % função: função a ser avaliada
              % xk    : Ponto onde se deseja avaliar a função
              [faprx2,dfobjnew] = feval(funcao,xksup);
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
                % Funcao gera aproximacoes para as funcoes para o metodo da Seção Aurea.
                [faprx1] = otapp_a2(Hess,XK,k,xkinf,fnew,dfnew);
            else
                % ==========================================
                % ==========================================
                % -------------------------------------------------
                % fob   : valor da função objetivo para a iteração.
                % dfobj : gradiente da função objetivo.
                % função: função a ser avaliada
                % xk    : Ponto onde se deseja avaliar a função
                [faprx1,dfobjnew] = feval(funcao,xkinf);
                 icfunc = icfunc + 1;
                % ==========================================
                % ==========================================
            end
           
          % Flag para impedir numero infinito de loops
            fim = fim + 1;
         end
  end
% ------------------------------------------                  
% Valor final para o alfa ótimo
  if faprx1 > faprx2
     limalfa = alfainf;
  else
     limbeta = alfasup;
  end
  alfaotim = abs(limalfa+limbeta)/2;
% ================================================================
% Gera aproximacao quadratica para a funcao.
  function [funaprox] = otapp_a2(Hess,XK,k,xkdavez,fnew,dfnew)
% Parâmetros de saída         
% fim     : Flag para a rotina.
% XK      % Historico dos xk's
% k       : Numero da iteração
% fnew    : Valor atual da função.
% dfnew   : Valor atual do gradiente da função.
% xkdavez : Valor atual de xk.
% -------------------------------------
% Gera uma aproximação para a funçao em xkdavez.
  dx = xkdavez - XK(:,k);
  funaprox = fnew + dfnew'*dx + 0.5*dx'*Hess*dx;
% =====================================================================