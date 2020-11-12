function dydx=centerdiff(x,y,x_res)
x1=[x(1)-x_res;x];
x2=[x;x(end)+x_res];
y1=interp1(x,y,x1,'linear','extrap');
y2=interp1(x,y,x2,'linear','extrap');
dydx=(diff(y1)+diff(y2))./(diff(x1)+diff(x2));
end