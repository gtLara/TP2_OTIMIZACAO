  clear;
  format compact
  pfiletype = '-dpng';
  pfileext  = '.png ';
% ========================================================
% ========================================================
  obj = @(x1,x2)     x1^2 + x2^2 - 2*x1 - 2*x2 + 2;
  
% % %   limx = 2;
% % %   limy = 3;
% % %   x1 = linspace(-limx,limx/2,100);
% % %   x2 = linspace(-limy,limy/3,100);
  limx = 4;
  limy = 4;
  x1 = linspace(0.01,4,91);
  x2 = linspace(0.01,4,95);

% Specified tolerance.
  epslon = 1.0e-08;
% -------------------------------------------
% Specified value considered zero = epslon^2;
  valzero = epslon^3;
% ------------------------------------
% Specified value considered infinite (eps = 2.2204e-016).
  INFMAX = 1/(eps^4);
  minfun  =  INFMAX;
  maxfun  = -INFMAX;
  mincon1 =  INFMAX;
  maxcon1 = -INFMAX;
  mincon2 =  INFMAX;
  maxcon2 = -INFMAX;
  mincon3 =  INFMAX;
  maxcon3 = -INFMAX;
  for i=1:length(x1)
      for j=1:length(x2)
        % =================================
          f(i,j)    = obj(x1(i),x2(j));
          if f(i,j) > maxfun
             maxfun = f(i,j);
          end
          if f(i,j) < minfun
             minfun = f(i,j);
          end
        % =================================
          con1(i,j) = g1(x1(i),x2(j));
          if con1(i,j) > maxcon1
             maxcon1 = con1(i,j);
          end
          if con1(i,j) < mincon1
             mincon1 = con1(i,j);
          end
        % =================================
          con2(i,j) = g2(x1(i),x2(j));
          if con2(i,j) > maxcon2
             maxcon2 = con1(i,j);
          end
          if con2(i,j) < mincon2
             mincon2 = con1(i,j);
          end  
        % =================================
          con3(i,j) = g3(x1(i),x2(j));
          if con3(i,j) > maxcon3
             maxcon3 = con3(i,j);
          end
          if con3(i,j) < mincon3
             mincon3 = con3(i,j);
          end         
      end
  end
  
% vetval = linspace(minfun,maxfun,20);
  delval = maxfun - minfun;
  vetval = [];
% n1 = 4;
  n1 = 3;
  for i1 = 1:n1
      divxy = 10^(n1-i1+1);
      for j1 = 1:n1
          vetval = [vetval, minfun + delval*(2^j1)/divxy] ;
      end
  end
  [c,h] = contour(x1,x2',f',vetval);
  xlabel('x_1')
  ylabel('x_2')
% clabel(c,h)
  axis equal 
 
% % %   limx = 2;
% % %   limy = 3;
% % %   x1 = linspace(-limx,limx/2,100);
% % %   x2 = linspace(-limy,limy/3,100);
  axis([0 4 0 4])
  c1 = ocontourc(x1,x2',con1',[0  1e6]);
  c2 = ocontourc(x1,x2',con2',[0  1e6]);
  c3 = ocontourc(x1,x2',con3',[0 -1e6]);
  
  hold on
  hatchedcontours(c1,'k');
  hatchedcontours(c2,'r');
  hatchedcontours(c3,'b');
  hold off
  print(pfiletype,'-r600',strcat('HatchedExample', pfileext));

