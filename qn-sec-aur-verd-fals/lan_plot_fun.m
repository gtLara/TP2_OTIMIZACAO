% ======================================================
% ======================================================
  function lan_plot_fun
  clear;
   
% ----------------------
% Função a ser estudada.
  icaso = 1;
  if icaso == 1
     funcao = 'fexlivro';
  end 
% ----------------------
% Função a ser estudada.
  if icaso == 2
     funcao = 'fun_rosensuzuki_irr';
  end
% --------------------------------------------

% ----------------------
% Função a ser estudada.
  if icaso == 3
     funcao = 'f_tiltednormcond';
  end

% -------------------------------------------------------------
% Definicao da dimensao do problema.
  disp('    '),
  disp('  --------------------------------------------------');
  disp('     dimensao do problema :  nxk   =  2 .'),
  nxk  =  2;
  disp('  -------------------------------------------------');
  disp('    '),
% -------------------------------------------------------------
% Vetores limites inferiores e superiores de x (Xmax and Xmin).
% ------------------- 
% Conjunto de dados
  x1 = linspace(-10,10,100);
  x2 = linspace(-10,10,100);
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
% saidas graficas
  figure(1), surf(X11,X22,X33); title('function obj');
  figure(2), surf(X11,X22,X44); title('gradient function objective');  
  end
% ======================================================
% ======================================================  