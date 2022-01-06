function [fobj,dfobj]=fteste3(x)
  
a=0.0263;

% Função objetivo.
dim = length(x);
A1 = diag(1:dim);
c1 = ones(dim,1);
fobj = 0.5*(x - c1)'*A1*(x - c1);

% Gradiente da função objetivo.
dfobj = A1*(x - c1); 

% % Corrige termos não quadráticos
% if a~=0
%    for i=1:dim   
%        fobj = fobj + a*(x(i,1)^3);  
%        dfobj(i,1) = dfobj(i,1) + 3*a*x(i,1)^2;      
%    end    
% end

% Corrige termos não quadráticos
if a~=0
   for i=1:dim   
       fobj = fobj + a*((x(i,1)-1)^3);  
       dfobj(i,1) = dfobj(i,1) + 3*a*(x(i,1)-1)^2;      
   end    
end