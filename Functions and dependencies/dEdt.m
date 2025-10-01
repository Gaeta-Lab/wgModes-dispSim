function nl_product=dEdt(Etemp,gamma_nl,L)

nl_product=1i*gamma_nl*abs(Etemp).^2*L.*Etemp;

% N=length(excited_modes);
% product=zeros(N,length(temp(1,:)));
%   for p=1:N
%       for q=1:N
%           for r=1:N
%               for s=1:N
%                   product(p,:)=product(p,:)+1i*gam_mat(p,q,r,s)*temp(q,:).*temp(r,:).*conj(temp(s,:));
%               end
%           end
%       end
%   end