file = fopen("experimento_2.txt", "w")
format = "%s %s %s %s %s\n";

% fexlivro

xvec = [[9;9], [-3;2]];
avec = [0.0263, -0.0263, 1]

for a=avec

    for aqn=0:3 % row grid
        [time, nit, nav, errx, errsol] = get_row(aqn, 2, xvec(:,1), a, 0);
        fprintf(file, format, time, nit, nav, errx, errsol);
    end

    fprintf(file, "\n");

    for aqn=0:3 % row grid
        [time, nit, nav, errx, errsol] = get_row(aqn, 2, xvec(:,2), a, 1);
        fprintf(file, format, time, nit, nav, errx, errsol);
    end

    fprintf(file, "\n\n");

end

fclose(file)
