function [p,d]=cheb_matrix_set(n)
%
nm=n-1;
for i=1:n 
    p(i,1)=cos(pi*(i-1)/nm);
end
%
d=zeros(n);
%
d(1,1)=(1.0+2.0*nm*nm)/6.0;
d(1,n)=(-1.0)^nm/(p(1)-p(n));
d(n,1)=(-1.0)^nm/(p(n)-p(1));
d(n,n)=-(1.0+2.0*nm*nm)/6.0;
%
for i=2:nm
    im=i-1;
    d(i,1)=(-1.0)^im/(2.0*(p(i)-p(1)));
    d(i,n)=(-1.0)^(im+nm)/(2.0*(p(i)-p(n)));
end
%
for j=2:nm
    jm=j-1;
    d(1,j)=(-1.0)^jm*2.0/(p(1)-p(j));
    d(n,j)=(-1.0)^(jm+nm)*2.0/(p(n)-p(j));
end  
%
for i=2:nm
    for j=2:nm
        if i==j
            d(i,i)=-p(i)/(2.0*(1.0-p(i)^2));
        else
            d(i,j)=(-1.0)^(i+j-2)/(p(i)-p(j));
        end
    end
end    