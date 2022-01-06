% Oh Senhor Jesus!!!

clc
clear

global a

dados = 50;
for aval=1:dados
    for pa=1:4
        % ALGORITMO QUASE NEWTON:
        % 0 --> Algoritmo DFP (tetha = 0.0 , pho = 1.0)
        % 1 --> Algoritmo BFGS(tetha = 1.0 , pho = 1.0)
        % 2 --> Adaptação de Huang
        % 3 --> Adaptação de Biggs
        Imetqn = pa-1;
        
        % PRECISÃO:
        epslon = 1e-6;
        
        % DEFINIÇÃO DA FUNÇÃO A SER ESTUDADA:
        funcao = 'fexlivro_ex6';
        
        % DIMENSÃO DO ESPAÇO DE BUSCA:
        dim = 2;
        
        % --------------------------------------------------------------------%
        % Limite superior do intervalo de busca
        x_sup = [10; 10];
        a_sup = max(roots([36*x_sup(1,1)*x_sup(2,1) 144*x_sup(2,1)+48*x_sup(1,1) 48]));
        % Limite inferior do intervalo de busca
        x_inf = [-10; -10];
        a_inf = min(roots([36*x_inf(1,1)*x_inf(2,1) 144*x_inf(2,1)+48*x_inf(1,1) 48]));
        % --------------------------------------------------------------------%
        
        % OPÇÕES PARA O VALOR DO PARÂMETRO a:
        a = 0.0263;
        
        %.....................................................................%
        %     Val_a = 0;                                                      %
        %     if Val_a == 0                                                   %
        %         % Definição do valor a = 0                                  %
        %         a = (a_sup + a_inf)/2;                                      %
        %     end                                                             %
        %     if Val_a == 1                                                   %
        %         % Definição do valor a positivo próximo de zero.            %
        %         a = (a_sup + 0.95*a_inf)/2;                                 %
        %     end                                                             %
        %     if Val_a == -1                                                  %
        %         % Definição do valor a negativo próximo de zero.            %
        %         a = (0.95*a_sup + a_inf)/2;                                 %
        %     end                                                             %
        %.....................................................................%
        
        % SOLUÇÃO ÓTIMA xstar:
        xstar(1,1) = min(roots([0.1875*a^3 3*a^2 (0.25*a^2+14*a) (2*a+4) ((1/12)*a+8/6)])); % -1/3;
        xstar(2,1) = min(roots([0.1875*a^3 a^2 ((4/3)*a+6*a) 4 2])); % -1/2
        
        % FUNÇÃO OBJETIVO DA SOLUÇÃO ÓTIMA xstar:
        fobjstar = fexlivro(xstar);
        
        % --------------------------------------------------------------------%
        % Expoente para definir os limites inferiores e superiores de x
        dmax = 1;
        % --------------------------------------------------------------------%
        
        % DEFINIÇÃO DOS LIMITES INF E SUP DE x (Xmax E xMin) E DO PONTO INICIAL x0:
        x0 = -10 + rand(dim,1)*(10-(-10));
        Xmin = [-10;-10];
        Xmax = [10;10];
        
        % DEFINIÇÃO DO NÚMERO MÁXIMO DE ITERAÇÕES:
        MAXITER = 2*dim*100;
        
        % DEFINIÇÃO DA OPERAÇÃO A SER REALIZADA DA SEÇÃO ÁUREA:
        % 0 --> Falsa seção aurea
        % 1 --> Verdadeira seção aurea
        isa_FV = 0;
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
            [XK1,Hfobj,k,kgold,icfunc] = otqnmat_a77(funcao,isa_FV,Imetqn,MAXITER,x0,Xmin,Xmax,epslon);
            tempo(exe) = toc;
            clc
        end
        tmedio(aval,pa) = mean(tempo);
        ni(aval,pa) = length(Hfobj - 1);
        na(aval,pa) = icfunc;
        err_per_x(aval,pa) = 100*norm(xstar - XK1(:,end))/max(1,norm(xstar));
        err_per_fobj(aval,pa) = 100*norm(fobjstar - Hfobj(end))/max(1,norm(fobjstar));
    end
end

% TEMPO MÉDIO:
media_t = mean(tmedio);
mediana_t = median(tmedio);
desvio_t = std(tmedio);
erro_t = std(tmedio)/sqrt(dados);
variancia_t = var(tmedio);
minimo_t = min(tmedio);
maximo_t = max(tmedio);
figure(1)
clf
boxplot(tmedio,{'DFP','BFGS','Huang','Biggs'})
xlabel('Métodos Quase-Newton')
ylabel('Tempo de processamento médio (s)')

% NÚMERO DE ITERAÇÕES:
media_ni = mean(ni);
mediana_ni = median(ni);
desvio_ni = std(ni);
erro_ni = std(ni)/sqrt(dados);
variancia_ni = var(ni);
minimo_ni = min(ni);
maximo_ni = max(ni);
figure(2)
clf
boxplot(ni,{'DFP','BFGS','Huang','Biggs'})
xlabel('Métodos Quase-Newton')
ylabel('Nº de Iterações')

% NÚMERO DE AVALIAÇÕES DE f(x):
media_na = mean(na);
mediana_na = median(na);
desvio_na = std(na);
erro_na = std(na)/sqrt(dados);
variancia_na = var(na);
minimo_na = min(na);
maximo_na = max(na);
figure(3)
clf
boxplot(na,{'DFP','BFGS','Huang','Biggs'})
xlabel('Métodos Quase-Newton')
ylabel('Nº de Avaliações de f(x)')

exe = 1:dados;
% ERRO DA SOLUÇÃO:
media_x = mean(err_per_x);
mediana_x = median(err_per_x);
desvio_x = std(err_per_x);
erro_x = std(err_per_x)/sqrt(dados);
variancia_x = var(err_per_x);
minimo_x = min(err_per_x);
maximo_x = max(err_per_x);
figure(4)
clf
plot(exe,err_per_x(:,1),'-*')
hold on
plot(exe,err_per_x(:,2),'-o')
plot(exe,err_per_x(:,3),'-d')
plot(exe,err_per_x(:,4),'-s')
legend('DFP','BFGS','Huang','Biggs')
xlabel('Execução')
ylabel('Erro da Solução (%)')

% ERRO DE f(x):
media_fobj = mean(err_per_fobj);
mediana_fobj = median(err_per_fobj);
desvio_fobj = std(err_per_fobj);
erro_fobj = std(err_per_fobj)/sqrt(dados);
variancia_fobj = var(err_per_fobj);
minimo_fobj = min(err_per_fobj);
maximo_fobj = max(err_per_fobj);
figure(5)
clf
plot(exe,err_per_fobj(:,1),'-*')
hold on
plot(exe,err_per_fobj(:,2),'-o')
plot(exe,err_per_fobj(:,3),'-d')
plot(exe,err_per_fobj(:,4),'-s')
legend('DFP','BFGS','Huang','Biggs')
xlabel('Execução')
ylabel('Erro da Função Objetivo (%)')

disp('TEMPO DE PROCESSAMENTO:')
disp(['média = ',num2str(media_t)])
disp(['mediana = ',num2str(mediana_t)])
disp(['desvio padrão = ',num2str(desvio_t)])
disp(['erro padrão = ',num2str(erro_t)])
disp(['variância = ',num2str(variancia_t)])
disp(['mínimo = ',num2str(minimo_t)])
disp(['máximo = ',num2str(maximo_t)])
disp('----------------------------------------------')
disp('Nº DE ITERAÇÕES:')
disp(['média = ',num2str(media_ni)])
disp(['mediana = ',num2str(mediana_ni)])
disp(['desvio padrão = ',num2str(desvio_ni)])
disp(['erro padrão = ',num2str(erro_ni)])
disp(['variância = ',num2str(variancia_ni)])
disp(['mínimo = ',num2str(minimo_ni)])
disp(['máximo = ',num2str(maximo_ni)])
disp('----------------------------------------------')
disp('Nº DE AVALIAÇÕES DE f(x):')
disp(['média = ',num2str(media_na)])
disp(['mediana = ',num2str(mediana_na)])
disp(['desvio padrão = ',num2str(desvio_na)])
disp(['erro padrão = ',num2str(erro_na)])
disp(['variância = ',num2str(variancia_na)])
disp(['mínimo = ',num2str(minimo_na)])
disp(['máximo = ',num2str(maximo_na)])
disp('----------------------------------------------')
disp('ERRO PERCENTUAL DE X:')
disp(['média = ',num2str(media_x)])
disp(['mediana = ',num2str(mediana_x)])
disp(['desvio padrão = ',num2str(desvio_x)])
disp(['erro padrão = ',num2str(erro_x)])
disp(['variância = ',num2str(variancia_x)])
disp(['mínimo = ',num2str(minimo_x)])
disp(['máximo = ',num2str(maximo_x)])
disp('----------------------------------------------')
disp('ERRO PERCENTUAL DE f(x):')
disp(['média = ',num2str(media_fobj)])
disp(['mediana = ',num2str(mediana_fobj)])
disp(['desvio padrão = ',num2str(desvio_fobj)])
disp(['erro padrão = ',num2str(erro_fobj)])
disp(['variância = ',num2str(variancia_fobj)])
disp(['mínimo = ',num2str(minimo_fobj)])
disp(['máximo = ',num2str(maximo_fobj)])