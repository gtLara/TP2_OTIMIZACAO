file = fopen("experimento_5.txt", "w")
format = "%s %s %s %s %s\n";

% fex1

x = [9;9]
for modusoperandi=1:-1:0

    for aqn=0:3 % row grid
        [time, nit, nav, errx, errsol] = get_row(aqn, 5, x, 0.0263, modusoperandi);
        fprintf(file, format, time, nit, nav, errx, errsol);
    end
    fprintf(file, "\n");
end

fclose(file)
