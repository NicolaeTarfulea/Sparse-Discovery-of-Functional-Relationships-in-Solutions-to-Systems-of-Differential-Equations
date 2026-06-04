function dxdt=enzyme(t,X)
k1=0.1;
km1=0.2;
k2=0.3;
dxdt=[k1*X(3)*X(2)-(km1+k2)*X(1); -k1*X(3)*X(2)+(km1+k2)*X(1); -k1*X(3)*X(2)+km1*X(1); k2*X(1)];
end

