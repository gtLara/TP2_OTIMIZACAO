% ------------------------------------------------------------------
% Algoritmos Quase Netwton - BFGS
  function [XK,Hfobj,k,kgold,icfunc]= otqnmat_a22(funcao,isa_FV,Imetqn,MAXITER,x0,Xmin,Xmax,epslon)
% ------------------------------------------------------------------
% For solving nonlinear optimization problems.
%     Minimize   fobj(x)
%	         xmin  <= x <=  xmax
% -----------------------------------------------------------------
% Iniciar variáveis.
  k = 1;             % Numero de iteracoes
  fim = 0;           % Flag para termino da execucao do algoritmo
  XK =[];            % Historico dos xk's
  Hfobj = [];        % Historico dos valores da função objetivo
  Hgrad = [];        % Historico dos valores do gradiente 
  xkold =[];         % valor de xk para  a iteracao anterior
  gkold =[];         % valor do gradiente para a iteracao anterior
  nxk = length(x0);  % @ DAGV -Numero de variaveis
  Hess = eye(nxk);   % Aproximaçao para a Hessiana (Identidade).   
  kgold= 0;          % numero real de avaliacoes da funcao.
% isa_FV = 0 Falsa seção aurea
% isa_FV = 1 verdadeira seção aurea;
% -----------------------
% Armazena o valor de x0.
  xk = x0;
  XK = [XK xk];

% Fator modificador dos limites de tolerância epslon.
  keps = 1.0e-03;

  
% Define o numero de avaliações reais da função objetivo.
  icfunc = 0;
  
% -----------------------
% Cria a matriz historico de dxk 
  Hdxk = zeros(nxk,3);
  
  while (~fim)
      % -------------------------------------------------
      % fob   : valor da função objetivo para a iteração.
      % dfobj : gradiente da função objetivo.
      % função: função a ser avaliada
      % xk    : Ponto onde se deseja avaliar a função
        icfunc = icfunc + 1;
        [fobjnew,dfobjnew] = feval(funcao,xk);
      % ------------------------------------
      % Armazena valor função objetivo.
        Hfobj = [Hfobj ; fobjnew];
      
      % ----------------------------------
      % Armazena o gradiente
        Hgrad=[Hgrad; dfobjnew'];
      
      % ----------------
      % Valor de deltaf.
        if k > 1
         % Verifica se a variação da função objetivo nula.
           if (abs(Hfobj(k) - Hfobj(k-1))) < epslon
              disp('     ')
              disp('     Valor da variação da função objetivo nulo.'),
              disp('     ')
              return
           end
        end
      
      % --------------------------------------------------
      % Norma vetor gradiente para função objetivo é nula.
        if norm(dfobjnew) <= keps*epslon
           disp('    Norma do gradiente função objetivo nula.  '),
           return
        end
      
      % ----------------------------
      % Define o valor de dxk em xk.
        dxk = -Imatqn*dfobjnew; % @ DAGV - Direcao de busca        
      
      % -------------------------------------------
      % Define direção de busca com vetor unitário
        mdxk = max(abs(dxk));
        if mdxk <= keps*epslon
           disp('    O módulo do vetor direção de busca nulo. '),
           return
        else
           dxk = dxk/mdxk;% @ dxk - Normalizacao.
        end
      
      % ----------------------------------------------------------
      % Valor do alfa ótimo para a função de mérito (Secao Aurea).
        if k > 1
           [Hess,IER]=otihmmv2_a2(tetha,pho,fobjold,fobjnew,dfobjold,...
                                  dfobjnew,xkold,xknew,Hess,1,corr);
           [alfaotim,icfunc] = otgoldsc_a22(funcao,icfunc,isa_FV,Hess,XK,...
            dxk,k,fobjold,fobjnew,dfobjold,dfobjnew,Xmin,Xmax,MAXITER);
        else
           if abs(dfobjnew'*dfobjnew) > 1.0e-3*epslon
              alfaotim = abs(0.05*fobjnew/abs(dfobjnew'*dfobjnew));
           else
              alfaotim = 0.025;
           end
        end
       
      % ---------------- ###############   ---------------------------------
      % ---------------- ###############   ---------------------------------
      % Gera a matriz aproximação para a hessiana (BFGS).
      % Determinar xknew, xkold, gknew e gkold.
        if k > 1 
         % Aproximação método BFGS.(tetha = 1.0; pho = 1.0).         
           corr=0;
           tetha = 1.0; % Se BFGS
           pho   = 1.0;    
           if Imetqn == 0 % Se DFP
              tetha = 0.0;
              pho = 1.0;
           elseif Imetqn == 2
              corr=1;
           elseif Imetqn==3
              corr=2;
           end         
           [Imatqn,IER]=otihmmv2_a2(tetha,pho,fobjold,fobjnew,dfobjold,...
                                    dfobjnew,xkold,xknew,Imatqn,-1,corr);
        else
         % Aproximaçao para a inversa da Hessiana(Identidade).
           Imatqn = eye(nxk);
        end
        
      % ----------------------------
      % Define o valor de dxk em xk.
        dxk = -Imatqn*dfobjnew; % @ DAGV - Direcao de busca 
      % ---------------- ###############   ---------------------------------
      % ---------------- ###############   ---------------------------------     
            
      % -------------------------------------------------
      % verifica se dxk esta oscilando. Se estiver reduz o alfa otimo 
      % por um fator de 10
        oscilacao = 0;
        for i=1:length(dxk)
            Hdxk(i,1) = Hdxk(i,2);
            Hdxk(i,2) = Hdxk(i,3);
            Hdxk(i,3) = dxk(i);
        end          
        ss = sign(Hdxk);
        for i=1:length(dxk)
            if ss(i,3) ~= ss(i,2)
               if ss(i,3)==ss(i,1)
                  oscilacao = 1;
              end                
            end 
        end
        if oscilacao == 1
           fator = 10*(1 - exp(-(k+1)/7));
           alfaotim = alfaotim/fator;
        end    
      
      % -----------------------------------------
      % Armazena xk em xkold e calcula o novo xk.
        xkold = xk;
        xk = xk + alfaotim*dxk;
      
      % ------------------------------------
      % Estuda viabilidade do novo vetor xk.
        for i =1:nxk
            if xk(i,1) < Xmin(i,1)
               xk(i,1) = Xmin(i,1);
            end
            if xk(i,1) > Xmax(i,1)
               xk(i,1) = Xmax(i,1);
            end
        end
      
      % ---------------------
      % Armazena xk em xknew.
        xknew = xk;
      
      % ---------------------
      % Armazena valor de xk.
        XK = [XK xk];
      
      % -----------------------------
      % Valor de deltax =(xknew - xkold).
        if norm(xknew - xkold) <= epslon
           disp('     Valor delta x nulo.'),
           return
        end
      
      % ---------------------------
      % disp(' Novo valor de xk '),
        xk'
        figure(10);
        clf;
        erxk = XK - repmat(xk,1,k+1);
        plot(sqrt(sum(erxk.*erxk,1)))
        drawnow
        disp( '                ' );
      % -----------------------------
     
      % --------------------------------------------------
      % Estocar fobjnew em fobjold e dfobjnew em dfobjold 
        fobjold = fobjnew;
        dfobjold = dfobjnew;
      
      % -------------------------------------------------
      % Verificar parada pelo número máximo de iterações.
        if k >= MAXITER
           fim = 1;
           disp('   Parada  número máximo de iterações.'),
        end
      
      % -------------- 
      % Nova iteração.
        if fim == 0
           k = k + 1;      
        end
  end % end while
% -------------------------------------------------