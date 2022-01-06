% Oh Senhor Jesus!!!

clc

% CARREGAR VARIÁVEIS PARA OS GRÁFICOS:
load DFP
load BFGS
load Huang
load Biggs

% SELECIONAR FUNÇÃO:
funcao = 'f_tiltednormcond';

% GRÁFICOS:
% Para fazer zoom nas figuras
zoom on;
[ll,cc]= size(XK1);
if ll <= 2
    figure(1);
    clf;
    % Plota o vetor x
    if cc > 1
        for i=1:length(x0)
            subplot(length(x0),1,i);
            plot(XK1_DFP(i,:),'-*');
            hold on
            plot(XK1_BFGS(i,:),'-o');
            plot(XK1_Huang(i,:),'-d');
            plot(XK1_Biggs(i,:),'-s');
            xlabel('Iteração'); ylabel(strcat('x_',num2str(i)));
            legend('DFP','BFGS','Huang','Biggs')
        end
    end
end
% Plota a função objetivo
figure(2);
clf;
plot(Hfobj_DFP,'-*');
hold on
plot(Hfobj_BFGS,'-o');
plot(Hfobj_Huang,'-d');
plot(Hfobj_Biggs,'-s');
xlabel('Iteração'); ylabel('f_o_b_j');
legend('DFP','BFGS','Huang','Biggs')

% Plota curvas de nível e soluções
if ll==2
%     x1 = linspace(0,3,40); % fex3
%     x2 = linspace(0,3,40); % fex3
    x1 = linspace(-10,10,40); % Outras f(x)
    x2 = linspace(-10,10,40); % Outras f(x)
    [X11,X22] = meshgrid(x1,x2);
    X33 = 0*X11;
    X44 = 0*X11;
    for i=1:length(x1)
        for j=1:length(x2)
            x = [X11(i,j); X22(i,j)];
            [fobj,dfobj] = feval(funcao,x);
            X33(i,j) = fobj;
            X44(i,j) = norm(dfobj);
        end
    end
    figure(3);
    clf;
    contour(X11,X22,X33,30);
    hold on
    plot(XK1_DFP(1,:),XK1_DFP(2,:),'-ob')
    plot(XK1_DFP(1,size(XK1_DFP,2)),XK1_DFP(2,size(XK1_DFP,2)),'*r')
    plot(XK1_DFP(1,size(XK1_DFP,2)),XK1_DFP(2,size(XK1_DFP,2)),'or')
    xlabel('x_1'); ylabel('x_2');
    figure(4);
    clf;
    contour(X11,X22,X33,30);
    hold on
    plot(XK1_BFGS(1,:),XK1_BFGS(2,:),'-ob')
    plot(XK1_BFGS(1,size(XK1_BFGS,2)),XK1_BFGS(2,size(XK1_BFGS,2)),'*r')
    plot(XK1_BFGS(1,size(XK1_BFGS,2)),XK1_BFGS(2,size(XK1_BFGS,2)),'or')
    xlabel('x_1'); ylabel('x_2');
    figure(5);
    clf;
    contour(X11,X22,X33,30);
    hold on
    plot(XK1_Huang(1,:),XK1_Huang(2,:),'-ob')
    plot(XK1_Huang(1,size(XK1_Huang,2)),XK1_Huang(2,size(XK1_Huang,2)),'*r')
    plot(XK1_Huang(1,size(XK1_Huang,2)),XK1_Huang(2,size(XK1_Huang,2)),'or')
    xlabel('x_1'); ylabel('x_2');
    figure(6);
    clf;
    contour(X11,X22,X33,30);
    hold on
    plot(XK1_Biggs(1,:),XK1_Biggs(2,:),'-ob')
    plot(XK1_Biggs(1,size(XK1_Biggs,2)),XK1_Biggs(2,size(XK1_Biggs,2)),'*r')
    plot(XK1_Biggs(1,size(XK1_Biggs,2)),XK1_Biggs(2,size(XK1_Biggs,2)),'or')
    xlabel('x_1'); ylabel('x_2');
end