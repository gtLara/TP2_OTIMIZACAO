file = fopen("experimento_4.txt", "w")
format = "%s %s %s %s %s\n";

% fex1

x = [-6;-5]
for modusoperandi=1:-1:0

    for aqn=0:3 % row grid
        [time, nit, nav, errx, errsol] = get_row(aqn, 4, x, 0.0263, modusoperandi);
        fprintf(file, format, time, nit, nav, errx, errsol);
    end
    fprintf(file, "\n");
end

fprintf(file, "\n\n");

x = - 10 + 20*rand(30, 1);
for modusoperandi=1:-1:0

    for aqn=0:3 % row grid
        [time, nit, nav, errx, errsol] = get_row(aqn, 1, x, 0.0263, modusoperandi);
        fprintf(file, format, time, nit, nav, errx, errsol);
    end
    fprintf(file, "\n");
end

fclose(file)
