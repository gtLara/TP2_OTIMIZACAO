  function [f,df] = f_tiltednormcond(x)
% printthis:
% Tilted and conditioned norm function:
% f = w ||A*x||_p + (w-1)<e1,A*x>
% parameters: A, w, p
% show: n, y, y
% pars.fgname = 'tiltednormcond'; % name of objective function (string)
% pars.p      = 2;                % parameter parsed to obj. function
% pars.w      = 8;                % parameter parsed to obj. function
% pars.nvar   = 10;               % number of variables
  p  = 2;
% w  = 8;
  w  = 2;

% % -----------------------------------
% % Definicao da operação a ser realizada.
%   disp('    '),
%   disp('  ------------------------------------------');
%   disp('    Number of variables:  pars.nvar.'),
%   pars.nvar = input('     Number of variables   = ');
%   disp('  -----------------------------------------');
%   disp('    '),
% % -------------------------------------------------------------

% pars.A      = rand(pars.nvar);  % parameter parsed to obj. function
% Hilbert matrix.
% HILB(N) is the N by N matrix with elements 1/(i+j-1),
% which is a famous example of a badly conditioned matrix.
  ndim = length(x);
  A  =  hilb(ndim);
% 
  Ax  = A*x;
  nAx = norm(Ax,p);
  f = w*nAx + (w-1)*Ax(1);
  df = w*x'*(A')*A /nAx + (w-1)*A(1,:);
  df = df';

