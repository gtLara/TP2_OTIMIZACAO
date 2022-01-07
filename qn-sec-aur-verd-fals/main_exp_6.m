pkg load statistics
file = fopen("experimento_6.txt", "w")
format = "%s %s %s %s\n";

% 0:direta  1:quadradica
opt = [1 0 1 0]; 
avec = [-0.0263, 0, 0, 0.0263];
N = 50;

for table=1:4
  
  % opt string (for plot title)
  if (opt(table) == 0)
    opt_str = "Aproximação Quadrática";
  else
    opt_str = "Avaliação Direta";
  endif

  % fill data matrices
  for i=1:N
    x= 20*rand(2,1) - 10;
    for aqn=0:3 % row grid
      % fexlivro
      [time, nit, nav, errx, errsol] = get_row(aqn, 2, x, avec(table), opt(table));
      time_vec(aqn+1,i) = str2num(time);
      nit_vec(aqn+1,i) = str2num(nit);
      nav_vec(aqn+1,i) = str2num(nav);
      errx_vec(aqn+1,i) = str2num(errx);
      errsol_vec(aqn+1,i) = str2num(errsol);    
    end
  end 

  % print data matrices to file in table format and generating graphs
  for data_id=1:5
    
    switch data_id
      case 1
        data = time_vec;
        title_str = "Tempo de Processamento Médio (s)";
      case 2
        data = nit_vec;
        title_str = "Número de Iterações";
      case 3
        data = nav_vec;
        title_str = "Número de Avaliações";
      case 4
        data = errx_vec;
        title_str = "Erro da solução";
      case 5
        data = errsol_vec;
        title_str = "Erro da função objetivo";
    endswitch
    
    % only generating tables and box plots for time, nit and nav
    if data_id <= 3
      for op=1:8
        row = zeros(1,4);
        for aqn=1:4
          switch op
            case 1
              row(aqn) = mean(data(aqn,:));
            case 2
              row(aqn) = median(data(aqn,:));
            case 3
              row(aqn) = std(data(aqn,:));
            case 4
              row(aqn) = std(data(aqn,:))/sqrt(length(data(aqn,:)));   
            case 5
              row(aqn) = var(data(aqn,:));
            case 6
              row(aqn) = min(data(aqn,:));
            case 7
              row(aqn) = max(data(aqn,:));
            case 8
              row(aqn) = length(data(aqn,:));         
          endswitch
          
        end
        fprintf(file, format, num2str(row(1)), num2str(row(2)), num2str(row(3)), num2str(row(4)));
      end
      fprintf(file, "\n");
      
      % generate box plots
      for aqn=1:4       
        
        if (aqn == 1)
          aqn_str = "DFP";
        elseif (aqn == 2)
          aqn_str = "BFGS";
        elseif (aqn == 3)
          aqn_str = "Huang";
        elseif (aqn == 3)
          aqn_str = "Biggs";
        end
        
        boxplot(data(aqn,:));
        title(sprintf("%s, %s, a = %s, %s", title_str, aqn_str, num2str(avec(table)), opt_str),'FontSize', 10);
        saveas(gcf, sprintf("plots/plot_%d_%d_%d.png", table, data_id, aqn));
      end
    else
      % generate box plots
      for aqn=1:4        
        if (aqn == 1)
          aqn_str = "DFP";
        elseif (aqn == 2)
          aqn_str = "BFGS";
        elseif (aqn == 3)
          aqn_str = "Huang";
        elseif (aqn == 3)
          aqn_str = "Biggs";
        end
        plot(data(aqn,:),'-o',2);
        title(sprintf("%s, %s, a = %s, %s", title_str, aqn_str, num2str(avec(table)), opt_str),'FontSize', 10);
        saveas(gcf, sprintf("plots/plot_%d_%d_%d.png", table, data_id, aqn));
      end
    endif
  end

  fprintf(file, "\n\n\n\n");
end

fclose(file)