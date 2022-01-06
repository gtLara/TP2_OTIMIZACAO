file = fopen("experimento_3.txt", "w")
format = "%s %s %s %s %s\n";

% fex1

dims = 5:5:25
for dim=dims
    x = ones(dim, 1) * 8;
    for aqn=0:3 % row grid
        [time, nit, nav, errx, errsol] = get_row(aqn, 1, x, 0.0263, 1);
        fprintf(file, format, time, nit, nav, errx, errsol);
    end
    fprintf(file, "\n");
end

fclose(file)
