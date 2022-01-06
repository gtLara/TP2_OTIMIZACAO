% Oh Senhor Jesus!!!

clc

x = linspace(-20,20,100);
a = 0;
for i=1:100
    fobj1(i) = -4*x(i)^2 + + 2*a*x(i)^3;
end
a = 0.0263;
for i=1:100
    fobj2(i) = -4*x(i)^2 + + 2*a*x(i)^3;
end
a = 0.0500;
for i=1:100
    fobj3(i) = -4*x(i)^2 + + 2*a*x(i)^3;
end
a = 0.2;
for i=1:100
    fobj4(i) = -4*x(i)^2 + + 2*a*x(i)^3;
end
a = 0.3;
for i=1:100
    fobj5(i) = -4*x(i)^2 + + 2*a*x(i)^3;
end

plot(x,fobj1,'--')
hold on
plot(x,fobj2,'.')
plot(x,fobj3,'-')
plot(x,fobj4,'*')
plot(x,fobj5,'o')
legend('a=0.0000','a=0.0263','a=0.0500','0=0.2000','a=0.3000')
xlabel('x')
ylabel('f_o_b_j')



