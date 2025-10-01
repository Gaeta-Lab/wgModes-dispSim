%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DO NOT CHANGE THIS FILE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dydx] = finite_diff_deriv(y,dx)

%COMPUTE A NUMERICAL DERIVATIVE USING 6TH ORDER CENTRAL DIFFERENCE FOR
%INTERIOR POINTS AND 6TH ORDER SINGLE-SIDED NEAR ENDPOINTS (ASSUMING EQUAL
%X SPACING OF dx).
%SOURCE: http://web.media.mit.edu/~crtaylor/calculator.html
n=length(y);
dydx=nan(n,1);
CD=[-1/60 3/20 -3/4 0 3/4 -3/20 1/60]/dx; %CENTRAL DIFFERENCE COEFFICIENTS
FD0=[-49/20 6 -15/2 20/3 -15/4 6/5 -1/6]/dx; %FORWARD DIFFERENCE
FD1=[-10 -77 150 -100 50 -15 2]/60/dx; %SEMI FORWARD DIFFERENCE 1
FD2=[2 -24 -35 80 -30  8 -1]/60/dx; %SEMI FORWARD DIFFERENCE 2
%INTERIOR
for k=4:n-3
    dydx(k)=dot(CD,y(k-3:k+3));
end
%NEAR ENDPOINTS
for k=1:3
    eval(['dydx(k)=dot(FD',num2str(k-1),',y(1:7));'])
    eval(['dydx(n-k+1)=dot(-flip(FD',num2str(k-1),'),y(n-6:n));'])
end
end