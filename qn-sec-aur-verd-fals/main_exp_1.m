file = fopen("experimento_1.txt", "w")
format = "%s %s %s %s %s\n";

% fex1
xvec = [[9;9], [-3;2], [-8;-6], [5;-7]];
i = 0;

for x=xvec
    for aqn=0:3 % row grid
        [time, nit, nav, errx, errsol] = get_row(aqn, 1, x, 0, mod(i, 2));
        fprintf(file, format, time, nit, nav, errx, errsol);
    end
    fprintf(file, "\n");
    i = i + 1;
end
fprintf(file, "\n\n");

% fexlivro

xvec = [[9;9], [-3;2]];
for aqn=0:3 % row grid
    [time, nit, nav, errx, errsol] = get_row(aqn, 2, xvec(:,1), 0, 0);
    fprintf(file, format, time, nit, nav, errx, errsol);
end
fprintf(file, "\n");
for aqn=0:3 % row grid
    [time, nit, nav, errx, errsol] = get_row(aqn, 2, xvec(:,2), 0, 1);
    fprintf(file, format, time, nit, nav, errx, errsol);
end
fprintf(file, "\n\n");

% fex3
xvec = [[9;9], [-3;2]];
for aqn=0:3 % row grid
    [time, nit, nav, errx, errsol] = get_row(aqn, 3, xvec(:,1), sqrt(2), 1);
    fprintf(file, format, time, nit, nav, errx, errsol);
end
fprintf(file, "\n");
for aqn=0:3 % row grid
    [time, nit, nav, errx, errsol] = get_row(aqn, 3, xvec(:,2), sqrt(2), 0);
    fprintf(file, format, time, nit, nav, errx, errsol);
end

fclose(file)
