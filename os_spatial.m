function output=os_spat_single(re,omega,n)
format short e

%Spatial Orr-Sommerfeld solver

%Actually a temporal Orr-Sommerfeld solver used within a secant shooting
%wrapper

%Temporal method: Supply real alpha, return complex omega
%Spatial method: Supply real omega, return complex alpha

%This method:
% 1. Supply real omega (omega_target)
% 2. Make crude guess at alpha that would return omega_target if input into
% temporal Orr-Sommerfeld equation (alpha1)
% 3. Solve temporal Orr-Sommerfeld equation for alpha1, returning omega1
% 4. Increment alpha1 slightly by small amount to provide alpha2, solve
% temporal Orr-Sommerfeld equation to return omega2
% 5. Now can compute d_omega/d_alpha and begin iterative method to
% determine omega_target. Final answer is the value of alpha that gives
% omega_target when supplied to the Orr-Sommerfeld solver

%Read boundary layer data 
bl_prof=bl_reader;

dimms=size(bl_prof);
nxp=dimms(1);
nyp=dimms(2)-1;

count=0;
istart=1;

%Initial alpha guess, to supply to solver
alpha=complex(omega/0.5,0);
omega=complex(omega,0);
alpha1=alpha;

%Target (supplied) omega
om_targ=omega;

        subplot(2,2,1)
           hold off
        subplot(2,2,2)
           hold off

           %Set discretisation
[p,d1]=cheb_matrix_set(n);  d2=d1^2;    d3=d2*d1;   d4=d3*d1;
    for k=1:nxp
    xloc(k)=bl_prof(k,nyp+1,1);
    end
    
         k=1;
           while xloc(k)<0.4
               k=k+1;
           end
         k=k+1;
         
%Set x-location here
xloc(k)         
    
    for j=1:nyp
        ul(j)=bl_prof(k,j,1);
        yl(j)=bl_prof(k,j,2);
    end

    p(1)=1.0-1.0e-6;
    L=0.1;
    
    fsj=int16(nyp*4/5);
    ufs=ul(fsj);

     for i=1:n
        y(i)=L*(1.0+p(i))./((1.0-p(i)));
     end
    
   %force freestream velocity to a chosen value using smooth function fy
        
   % fymax=8/log(1+exp(8))
   % fymin=log(1+exp(-8))/8*fymax
    for j=1:fsj
    fx(j)=(2*yl(j)-yl(fsj))/yl(fsj);
 %   fy(j)=(log(1+exp(8*fx(j)))/8);
    fy(j)=(exp(8*fx(j))/(1+exp(8*fx(j))));
    end
    
    fy=fy-fy(1);
    fy=fy/fy(fsj);
    
    for j=1:nyp
        if j<fsj
            ul(j)=ul(j)*(1-fy(j))+fy(j)*ufs;
        else
            ul(j)=ufs;
        end
        
    end
    
    subplot(2,2,4)
           hold on
           plot(fx,fy)
           hold on
           
    
    eta1=(-((y-L)./(L+y).^2)+1.0./(L+y));
    eta2=(2.0*(y-L)./(L+y).^3-2.0./(L+y).^2);
    eta3=(-6.0*(y-L)./(L+y).^4+6.0./(L+y).^3);
    eta4=(24.0*(y-L)./(L+y).^5-24.0./(L+y).^4);
    for i=1:n
        for j=1:n
            d4(i,j)=eta1(i)^4.0*d4(i,j)+6.0*eta1(i)^2*eta2(i)*d3(i,j)+(3.0*eta2(i)^2+4.0*eta1(i)*eta3(i))*d2(i,j)+eta4(i)*d1(i,j);
            d3(i,j)=eta1(i)^3*d3(i,j)+3.0*eta1(i)*eta2(i)*d2(i,j)+eta3(i)*d1(i,j);
            d2(i,j)=eta1(i)^2*d2(i,j)+eta2(i)*d1(i,j);
            d1(i,j)=eta1(i)*d1(i,j);
        end
    end    
    %jfree=int16(nyp*3/4)
    
    for i=1:n
        if y(i)<yl(fsj)
    %        f=(1-exp(-y(i)))*(1+cos((y(i)-yl(jfree))/yl(jfree)*pi()))*0.5;
    %        u(i)=spline(yl,ul,y(i))*(1-f)+f*ul(jfree);

        u(i)=spline(yl,ul,y(i));
  %       u(i)=(spline(yl,ul,y(i)))*(1-f)+f(
        else
           u(i)=ufs;
        end
    end

    % compute second derivative of base flow
    d2u=d2*u';

%COMPUTE OMEGA FOR FIRST ALPHA GUESS, ALPHA1
    
    diff=5
    alpha=alpha1
    % set Orr-Sommerfeld equation in matrix form A*vhat=omega*B*vhat
    ci=complex(0.0,1.0);

    for i=1:n
        for j=1:n
            A(i,j)=(u(i)-2.0*ci*alpha/re)*d2(i,j)+ci/(re*alpha)*d4(i,j);
            B(i,j)=d2(i,j);
            if i==j
                A(i,j)=A(i,j)+ci*alpha^3/re-u(i)*alpha^2-d2u(i);
                B(i,j)=B(i,j)-alpha^2;
            end
        end
    end

    for j=1:n
        A(1,j)=0.0;
        A(2,j)=0.0;
        A(n-1,j)=0.0;    
        A(n,j)=0.0;
        B(1,j)=0.0;
        B(2,j)=d1(1,j);
        B(n-1,j)=d1(n,j);
        B(n,j)=0.0;    
    end
    
    B(1,1)=1.0;
    B(n,n)=1.0;

    e=eig(inv(B)*A);
    eval=complex(0.0,-1.e10);
    for i=1:n
        if imag(e(i))>imag(eval)
            if real(e(i))<0.9
            eval=e(i);
            end
        end     
    end

    inc=complex(1e-3,1e-3)

    %COMPUTE OMEGA FOR INCREMENTED ALPHA, ALPHA2
    
    omega1=eval*alpha1;
    alpha2=alpha1+inc;
    alpha=alpha2;
    
        for i=1:n
        for j=1:n
            A(i,j)=(u(i)-2.0*ci*alpha/re)*d2(i,j)+ci/(re*alpha)*d4(i,j);
            B(i,j)=d2(i,j);
            if i==j
                A(i,j)=A(i,j)+ci*alpha^3/re-u(i)*alpha^2-d2u(i);
                B(i,j)=B(i,j)-alpha^2;
            end
        end
    end

    for j=1:n
        A(1,j)=0.0;
        A(2,j)=0.0;
        A(n-1,j)=0.0;    
        A(n,j)=0.0;
        B(1,j)=0.0;
        B(2,j)=d1(1,j);
        B(n-1,j)=d1(n,j);
        B(n,j)=0.0;    
    end
    
    B(1,1)=1.0;
    B(n,n)=1.0;

    e=eig(inv(B)*A);
    eval=complex(0.0,-1.e10);
    for i=1:n
        if imag(e(i))>imag(eval)
            if real(e(i))<0.9
            eval=e(i);
            end
        end     
    end

    omega2=eval*alpha2;
    
    iter=0;
    
    if omega2-omega1==0
        cph_r(k)=0
        alpha=0
    else

    while diff>1e-4
    iter=iter+1
    
        %COMPUTE D_ALPHA/D_OMEGA
        % START SEARCHING FOR OMEGA_TARGET
        
    dadw=(alpha2-alpha1)/(omega2-omega1);
    alpha=(om_targ-omega2)*dadw+alpha2;
    da=alpha-alpha2;
    
    for i=1:n
        for j=1:n
            A(i,j)=(u(i)-2.0*ci*alpha/re)*d2(i,j)+ci/(re*alpha)*d4(i,j);
            B(i,j)=d2(i,j);
            if i==j
                A(i,j)=A(i,j)+ci*alpha^3/re-u(i)*alpha^2-d2u(i);
                B(i,j)=B(i,j)-alpha^2;
            end
        end
    end

    for j=1:n
        A(1,j)=0.0;
        A(2,j)=0.0;
        A(n-1,j)=0.0;    
        A(n,j)=0.0;
        B(1,j)=0.0;
        B(2,j)=d1(1,j);
        B(n-1,j)=d1(n,j);
        B(n,j)=0.0;    
    end
    
    B(1,1)=1.0;
    B(n,n)=1.0;

    e=eig(inv(B)*A);
    eval=complex(0.0,-1.e10);
    for i=1:n
        if imag(e(i))>imag(eval)
            if real(e(i))<0.9
            eval=e(i);
            end
        end     
    end
    
    alpha1=alpha2;
    omega1=omega2;
    alpha2=alpha;
    omega2=alpha2*eval;
    diff=sqrt((imag(eval*alpha))^2);
    omega=omega2;
    
    if real(eval)<1e-8 
       diff=0
    elseif omega2-omega1==0
       diff=0       
    end
    
         end

    end 
    
    cph_r(k)=real(eval);
  %  cph_i(k)=imag(eval);

       if mod(k,10)==0
           for i=1:n/3
               u2(i)=d2u(n+1-i);
               y2(i)=y(n+1-i);
           end
           subplot(2,2,1)
           hold on
           plot(u2,y2)
           hold on
           
                  for i=1:n/2
               u2(i)=u(n+1-i);
               y2(i)=y(n+1-i);
           end
           subplot(2,2,2)
           hold on
           plot(u2,y2)
           hold on
           
       end

      %data is output variable when used as a subroutine
data(k,1)=xloc(k);

if cph_r(k)<1e-8 
   data(k,2)=0;
   data(k,3)=0;
else 
data(k,2)=(imag(alpha));
data(k,3)=cph_r(k);
end

if cph_r(k)>0
    if count==0
    istart=k;
    count=1;
    end
end

xloc(k)

istart=istart+1;
data(1,4)=istart;

output(2)=(alpha);
output(1)=(omega);

alpha=alpha
omega=omega

%subplot(1,1,1)
hold off

