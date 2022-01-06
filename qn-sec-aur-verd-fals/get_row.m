% Oh Senhor Jesus!!!

function [tmedio, ni, icfunc, err_per_x, err_per_fobj]=get_row(Imetqn, icaso, x0, a, isa_FV)

% ALGORITMO QUASE NEWTON:
% Imetqn = 0 --> Algoritmo DFP (tetha = 0.0 , pho = 1.0)
% Imetqn = 1 --> Algoritmo BFGS(tetha = 1.0 , pho = 1.0)
% Imetqn = 2 --> Adaptação de Huang
% Imetqn = 3 --> Adaptação de Biggs
% Imetqn = 0;
if Imetqn==0
    AQN = 'Algoritmo Quase-Newton: DFP';
elseif Imetqn==1
    AQN = 'Algoritmo Quase-Newton: BFGS';
elseif Imetqn==2
    AQN = 'Algoritmo Quase-Newton: Adaptação Huang';
elseif Imetqn==3
    AQN = 'Algoritmo Quase-Newton: Adaptação Biggs';
end    

% PRECISÃO:
epslon = 1e-6;

% DEFINIÇÃO DA FUNÇÃO A SER ESTUDADA:
% teste 1 --> fex1(a=0), fexlivro(a=0) e fex3
% teste 2 --> fexlivro
% teste 3 --> fex3
% teste 4 --> fun_rosensuzuki_irr e fex1 (a~=0)
% teste 5 --> f_tiltednormcond
% teste 6 --> fexlivro
%-------------------------------------------------
% icaso = 1 --> fex1
% icaso = 2 --> fexlivro
% icaso = 3 --> fex3
% icaso = 4 --> fun_rosensuzuki_irr
% icaso = 5 --> f_tiltednormcond
% icaso = 1;
if icaso==1
    F = 'Função objetivo: fex1';
elseif icaso==2
    F = 'Função objetivo: fexlivro';
elseif icaso==3
    F = 'Função objetivo: fex3';
elseif icaso==4
    F = 'Função objetivo: fun_rosensuzuki_irr';
elseif icaso==5
    F = 'Função objetivo: f_tiltednormcond';
end    

dim = size(x0)(1);

if icaso == 1
    funcao = 'fex1';
    
    % DIMENSÃO DO ESPAÇO DE BUSCA:
    
    %---------------------------------------------------------------------%
    % % Limite superior do intervalo de busca
    % x_sup = [10; 10];
    % a_sup = max(roots([36*x_sup(1,1)*x_sup(2,1) 144*x_sup(2,1)+48*x_sup(1,1) 48]));
    % % Limite inferior do intervalo de busca
    % x_inf = [-10; -10];
    % a_inf = min(roots([36*x_inf(1,1)*x_inf(2,1) 144*x_inf(2,1)+48*x_inf(1,1) 48]));
    %---------------------------------------------------------------------%
    
    % OPÇÕES PARA O VALOR DO PARÂMETRO a:
    % a = 0.0263; % a = 0 no enunciado
%.........................................................................%
%     Val_a = 0;                                                          %
%     if Val_a == 0                                                       %
%         % Definição do valor a = 0                                      %
%         a = (a_sup + a_inf)/2;                                          %
%     end                                                                 %
%     if Val_a == 1                                                       %
%         % Definição do valor a positivo próximo de zero.                %
%         a = (a_sup + 0.95*a_inf)/2;                                     %
%     end                                                                 %
%     if Val_a == -1                                                      %
%         % Definição do valor a negativo próximo de zero.                %
%         a = (0.95*a_sup + a_inf)/2;                                     %
%     end                                                                 %
%.........................................................................%

    % SOLUÇÃO ÓTIMA xstar:
%     xstar(1,1) = 1; % min(roots([0.1875*a^3 3*a^2 (0.25*a^2+14*a) (2*a+4) ((1/12)*a+8/6)]));
%     xstar(2,1) = 1; % min(roots([0.1875*a^3 a^2 ((4/3)*a+6*a) 4 2]));
    xstar = ones(dim,1);
    
    % FUNÇÃO OBJETIVO DA SOLUÇÃO ÓTIMA xstar:
    fobjstar = fex1(xstar, a);
    
    % --------------------------------------------------------------------%
    % Expoente para definir os limites inferiores e superiores de x
    dmax = 1;
    % --------------------------------------------------------------------%
end

if icaso == 2
    funcao = 'fexlivro';
    
    % DIMENSÃO DO ESPAÇO DE BUSCA:
    % dim = 2;
    
    % --------------------------------------------------------------------%
    % Limite superior do intervalo de busca
    % x_sup = [10; 10];
    % a_sup = max(roots([36*x_sup(1,1)*x_sup(2,1) 144*x_sup(2,1)+48*x_sup(1,1) 48]));    
    % % Limite inferior do intervalo de busca
    % x_inf = [-10; -10];
    % a_inf = min(roots([36*x_inf(1,1)*x_inf(2,1) 144*x_inf(2,1)+48*x_inf(1,1) 48]));
    % --------------------------------------------------------------------%
    
    % OPÇÕES PARA O VALOR DO PARÂMETRO a:
    % a = 1;
    
%.........................................................................%
%     Val_a = 0;                                                          %
%     if Val_a == 0                                                       %
%         % Definição do valor a = 0                                      %
%         a = (a_sup + a_inf)/2;                                          %
%     end                                                                 %
%     if Val_a == 1                                                       %
%         % Definição do valor a positivo próximo de zero.                %
%         a = (a_sup + 0.95*a_inf)/2;                                     %
%     end                                                                 %
%     if Val_a == -1                                                      %
%         % Definição do valor a negativo próximo de zero.                %
%         a = (0.95*a_sup + a_inf)/2;                                     %
%     end                                                                 %
%.........................................................................%
          
    % SOLUÇÃO ÓTIMA xstar:
    xstar(1,1) = min(roots([0.1875*a^3 3*a^2 (0.25*a^2+14*a) (2*a+4) ((1/12)*a+8/6)])); % -1/3; 
    xstar(2,1) = min(roots([0.1875*a^3 a^2 ((4/3)*a+6*a) 4 2])); % -1/2
    
    % FUNÇÃO OBJETIVO DA SOLUÇÃO ÓTIMA xstar:
    fobjstar = fexlivro(xstar, a);
    
    % --------------------------------------------------------------------%
    % Expoente para definir os limites inferiores e superiores de x
    dmax = 1;
    % --------------------------------------------------------------------%
end

if icaso == 3
    funcao = 'fex3';
    
    % DIMENSÃO DO ESPAÇO DE BUSCA:
    
    % --------------------------------------------------------------------%
    % Limite superior do intervalo de busca
    x_sup = [3; 3];        
    % Limite inferior do intervalo de busca
    x_inf = [0; 0];
    % --------------------------------------------------------------------%
    
    % SOLUÇÃO ÓTIMA xstar:
    xstar(1,1) = 1;
    xstar(2,1) = 1;
    
    % FUNÇÃO OBJETIVO DA SOLUÇÃO ÓTIMA xstar:
    fobjstar = fex3(xstar, a);
    
    % --------------------------------------------------------------------%
    % Expoente para definir os limites inferiores e superiores de x
    dmax = 1;
    % --------------------------------------------------------------------%
end

if icaso == 4
    funcao = 'fun_rosensuzuki_irr';
    
    % DIMENSÃO DO ESPAÇO DE BUSCA:
    
    % SOLUÇÃO ÓTIMA xstar:
    for i1 =1:dim
        xstar(i1,1) = (-1)^i1;
    end
    
    % FUNÇÃO OBJETIVO DA SOLUÇÃO ÓTIMA xstar:
    fobjstar = fun_rosensuzuki_irr(xstar);
    
    %---------------------------------------------------------------------%
    % Vetores limites inferiores e superiores de x (Xmax and Xmin)
    dmax = 4;
    %---------------------------------------------------------------------%
end

if icaso == 5
    funcao = 'f_tiltednormcond';
    
    % DIMENSÃO DO ESPAÇO DE BUSCA:
    
    % SOLUÇÃO ÓTIMA xstar:
    for i1 =1:dim
        xstar(i1,1) = 0.0;
    end
    
    % FUNÇÃO OBJETIVO DA SOLUÇÃO ÓTIMA xstar:
    fobjstar = f_tiltednormcond(xstar);
    
    %---------------------------------------------------------------------%
    % Vetores limites inferiores e superiores de x (Xmax and Xmin)
    dmax = 0.3010;
    %---------------------------------------------------------------------%
end

% DEFINIÇÃO DOS LIMITES INF E SUP DE x (Xmax E xMin) E DO PONTO INICIAL x0:
% x0 = [9;9];
Xmin = ones(size(x0)) * -10;
Xmax = ones(size(x0)) * 10;

%.........................................................................%
% for i1 =1:dim                                                           %
%     Xmin(i1,1) =-1*10.0^dmax;                                           %
%     Xmax(i1,1) = 1*10.0^dmax;                                           %
%     x0(i1,1) = Xmin(i1,1) + rand(1)*(Xmax(i1,1)-Xmin(i1,1));            %
% end                                                                     %
%.........................................................................%

% DEFINIÇÃO DO NÚMERO MÁXIMO DE ITERAÇÕES:
MAXITER = 2*dim*100;

% DEFINIÇÃO DA OPERAÇÃO A SER REALIZADA DA SEÇÃO ÁUREA:
% 0 --> Falsa seção aurea
% 1 --> Verdadeira seção aurea
if isa_FV==0
    TEC = 'Aproximação quadrática de f(x) em cada iteração';
else
    TEC = 'Técnica da seção áurea feita através da avaliação direta de f(x)';
end

% RESOLUÇÃO E RESULTADOS:
% funcao --> Função a ser aproximada
% Imetqn --> Método quase Newton Escolhido
% MAXITER --> Número máximo de iteracões
% x0 --> Ponto inicial
% [Xmin Xmax] --> Limites da função
% epslon --> Precisão
% Hfobj --> Armazena valor função objetivo
% XK1 --> Armazena as variáveis
for exe=1:10
    tic
    [XK1,Hfobj,k,kgold,icfunc] = otqnmat_a77(funcao,isa_FV,Imetqn,MAXITER,x0,Xmin,Xmax,epslon, a);
    tempo(exe) = toc;
    clc
end
tmedio = mean(tempo);
ni = length(Hfobj - 1);
err_per_x = 100*norm(xstar - XK1(:,end))/max(1,norm(xstar));
err_per_fobj = 100*norm(fobjstar - Hfobj(end))/max(1,norm(fobjstar));
% disp(AQN)
% disp(F)
% disp(TEC)
% disp(['Precisão = ',num2str(epslon)])
% % fprintf('Ponto inicial = %f\n',x0)
% disp(['Limite inferior = ',num2str(Xmin(1))])
% disp(['Limite superior = ',num2str(Xmax(1))])
% disp(['Número máximo de iterações = ',num2str(MAXITER)])
% disp(['Tempo de processamento = ',num2str(tmedio)])
% disp(['Número final de iterações = ',num2str(ni)])
% disp(['Número final de avaliações de fobj = ',num2str(icfunc)])
% disp(['Erro percetual fobj = ',num2str(err_per_fobj)])
% disp(['Erro percetual xsol = ',num2str(err_per_x)])

tmedio = num2str(tmedio)
ni = num2str(ni)
icfunc = num2str(icfunc)
err_per_fobj = num2str(err_per_fobj)
err_per_x = num2str(err_per_x)

% SALVAR VARIÁVEIS PARA GERAR GRÁFICOS:
% if Imetqn==0 % DFP
%     XK1_DFP = XK1;
%     Hfobj_DFP = Hfobj;
%     save DFP XK1_DFP Hfobj_DFP
% elseif Imetqn==1 % BFGS
%     XK1_BFGS = XK1;
%     Hfobj_BFGS = Hfobj;
%     save BFGS XK1_BFGS Hfobj_BFGS
% elseif Imetqn==2 % Huang
%     XK1_Huang = XK1;
%     Hfobj_Huang = Hfobj;
%     save Huang XK1_Huang Hfobj_Huang
% elseif Imetqn==3 % Biggs
%     XK1_Biggs = XK1;
%     Hfobj_Biggs = Hfobj;
%     save Biggs XK1_Biggs Hfobj_Biggs 
% end
