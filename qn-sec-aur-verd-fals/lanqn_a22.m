  clear;
% -----------------------------------------
% Aproximação Quadrática para a função objetivo.
% Imetqn : Controla a escolha do algoritmo DFP ou BFGS.
%      0 : Algoritmo DFP (tetha = 0.0 , pho = 1.0).
%      1 : Algoritmo BFGS(tetha = 1.0 , pho = 1.0).
%      2 : Huang.
%      3 : Biggs.
% -----------------------------------------
  disp ( '      ' );
  disp('  -------------------------------------------------');
  disp('    Controla a escolha do algoritmo DFP ou BFGS.'),
  Imetqn =input('     Algoritmo DFP [0]; Algoritmo BFGS [1]; Huang[2]; Biggs [3] =  ');
  disp('  -------------------------------------------------');
  disp ( '      ' );
% ---------------------------------------------------

% ------------------------------------
% Definição dos limites de tolerância.
% epslon = 1e-6; % gradiente
  epslon = 1e-13; % gradiente

% -------------------------------
% Define a função a ser estudada.
  icaso = 2;
% ----------------------
% Função a ser estudada.
  if icaso == 1
     funcao = 'fexlivro';
   % otimun=[-1/3; -1/2];
   % Dimension do espaco de busca.
     dim = 2;   
   % ===================================================
   % ===========================================================
   % Limite superior do intervalo de busca 
     x_sup = [10; 10]; 
     a_sup = max(roots([36*x_sup(1,1)*x_sup(2,1) 144*x_sup(2,1)+48*x_sup(1,1) 48]));
     
   % Limite inferior do intervalo de busca
     x_inf = [-10; -10]; 
     a_inf = min(roots([36*x_inf(1,1)*x_inf(2,1) 144*x_inf(2,1)+48*x_inf(1,1) 48]));
     Val_a = 0;
     if Val_a == 0
      % Definição do valor a = 0
        a = (a_sup + a_inf)/2;
     end
     if Val_a == 1
      % Definição do valor a positivo próximo de zero.
        a = (a_sup + 0.95*a_inf)/2;
     end
     if Val_a == -1
      % Definição do valor a negativo próximo de zero.
        a = (0.95*a_sup + a_inf)/2;
     end
   % xstar valor ótimo
     xstar(1,1) = min(roots([0.1875*a^3 3*a^2 (0.25*a^2+14*a) (2*a+4) ((1/12)*a+8/6)]));
     xstar(2,1) = min(roots([0.1875*a^3 a^2 ((4/3)*a+6*a) 4 2]));
   % ===========================================================
   % =================================================== 
   % optimun = [-1/3; -1/2];
   % xstar = [-1/3 ; -1/2];
   % -------------------------------------------------------------
   % Expoente para definir os limites inferiores e superiores de x.
     dmax = 1;
  end
  
% ----------------------
% Função a ser estudada.
  if icaso == 2
     funcao = 'fun_rosensuzuki_irr';
   % -----------------------------------
   % Definicao da operação a ser realizada.
     disp('    '),
     disp('  ------------------------------------------');
     disp('    Definição da dimensao de xk.'),
     dim = input('    dim    = ');
     disp('  -----------------------------------------');
     disp('    '),
   % -------------------------------------------------------------
   % otimun=[-1;1;-1;1;.....;-1;1];
     for i1 =1:dim
         xstar(i1,1) = (-1)^i1;
     end
   % Vetores limites inferiores e superiores de x (Xmax and Xmin).
     dmax = 4;
  end
   
% ----------------------
% Função a ser estudada.
  if icaso == 3
     funcao = 'f_tiltednormcond';
   % -----------------------------------
   % Definicao da operação a ser realizada.
     disp('    '),
     disp('  ------------------------------------------');
     disp('    Definição da dimensao de xk.'),
     dim = input('    dim    = ');
     disp('  -----------------------------------------');
     disp('    '),
   % -------------------------------------------------------------
   % otimun=[0;0;0;.....; 0];
     for i1 =1:dim
         xstar(i1,1) = 0.0;
     end

   % -------------------------------------------------------------
   % Expoente para definir os limites inferiores e superiores de x.
     dmax = 0.3010;   
  end

% -------------------------------------------------------------
% Vetores limites inferiores e superiores de x (Xmax and Xmin).
% O vetor x0 deve comecar em qualquer ponto do espaco de busca.
% ---------------------------
  for i1 =1:dim
      Xmin(i1,1) =-1*10.0^dmax;
      Xmax(i1,1) = 1*10.0^dmax;
      x0(i1,1) = Xmin(i1,1) + rand(1)*(Xmax(i1,1)-Xmin(i1,1));
  end

% ------------------------------------------
% Definição do número máximo de iterações.
  disp ( '      ' );
  disp('  ------------------------------------------  ');
  MAXITER = 2*dim*100;
  disp(['   Número máximo de iterações  =  ',num2str(MAXITER)]),
  disp('  -----------------------------------------   ');
  disp ( '      ' );
% ---------------------------------------------------
% isa_FV = 0 Falsa seção aurea
% isa_FV = 1 verdadeira seção aurea;
% -----------------------------------
% Definicao da operação a ser realizada.
  disp('    '),
  disp(' ----------------------------------------------------');
  disp('   isa_FV = 0  : Falsa seção aurea.        '),
  disp('   isa_FV = 1  : Verdadeira seção aurea.   '),
  isa_FV =input('    seção aurea Falsa ou Verdadeira    = ');
  disp('  ---------------------------------------------------');
  disp('    '),
% -------------------------------------------------------------

% -------------------------------------------------------
% funcao - funcao a ser aproximada
% Imetqn - metodo quase Newton Escolhido
% MAXITER - numero maximo de iteracoes
% x0 - ponto inicial
% [Xmin Xmax] - limites da funcao
% epslon - tolerancia
% Hfobj : Armazena valor função objetivo.
% XK1 : Armazena as variaveis
  [XK1,Hfobj,k,kgold,icfunc] = otqnmat_a77(funcao,isa_FV,Imetqn,MAXITER,x0,Xmin,Xmax,epslon);
% -------------------------------------------------------

% --------------------------------------------
  disp('  -------------------------------------------'); 
  disp(['   Número final de iterações  =  ',num2str(length(Hfobj)-1)]),
  disp('  -------------------------------------------');
  disp('    '),
% --------------------------------------------
% --------------------------------------------
  disp('  -------------------------------------------'); 
  disp(['   Número final de avaliações de fobj =  ',num2str(icfunc)]),
  disp('  -------------------------------------------');
  disp('    '),
% --------------------------------------------
% Avaliar o erro percentual xsol   
  err_per = 100*norm(xstar - XK1(:,end))/max(1,norm(xstar));
% --------------------------------------------
  disp('  -------------------------------------------'); 
  disp(['  erro percentual  xsol  =  ',num2str(err_per)]),
  disp('  -------------------------------------------');
  disp('    '),
% --------------------------------------------

% ----------------------------
% Para fazer zoom nas figuras.
  zoom on;
  [ll,cc]= size(XK1);
  if ll <= 2
     figure(7);
     clf;
   % Plota o vetor x.
     if cc > 1
        for i=1:length(x0)
            subplot(length(x0),1,i);
            plot(XK1(i,:)/abs(XK1(i,cc-1)));
        end  
     end
  end
% -----------------------------

% ------------------------
% Plota a função objetivo.
  figure(8);
  clf;
  plot(Hfobj/abs(Hfobj(cc-1)));
% -----------------------------