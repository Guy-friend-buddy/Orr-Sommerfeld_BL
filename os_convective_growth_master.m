function os_spatial_foil;
%This is the master file (i.e. the executable)
%This program will read in boundary layer data and then compute the spatial
%growth (i.e. as a function of x) of primary instabilities with
%user-specfiied frequencies

%Reynolds number
re=100000

%Number of modes (i.e. dimension of matrix) within Orr-Sommerfeld problem
nmodes=220

% if n==0
%     interval=0
%     nmax=1
% else
%     interval=(om2-om1)/(n);
%     nmax=n+1
% end

%Prepare plot windows
subplot (1,1,1)
hold off
subplot (2,2,1)
set(gca,'LineStyleOrder', '-*|:|o')
hold all
subplot (2,2,2)
set(gca,'LineStyleOrder', '-*|:|o')
hold all
subplot (2,2,3)

set(gca,'LineStyleOrder', '-*|:|o')
hold all
subplot (2,2,4)

set(gca,'LineStyleOrder', '-*|:|o')
hold all

%Read in boundary layer data via subroutine
bl_prof=bl_reader_subroutine;
dimms=size(bl_prof);
nxp=dimms(1);
nyp=dimms(2)-1;

%Hard-coded omega values
%Could also loop over an interval with equal spacing if preferred
nmax=4;
omstore(1)=24.38;
omstore(2)=21.24;
omstore(3)=26.7;
omstore(4)=24.25;

for i=1:nmax
%     Generate omega values automatically if preferred
%     om=om1+interval*(i-1);
%     omstore(i)=om;

    om=omstore(i);

    % Run Orr-Sommerfeld solver for specified omega
    % This subroutine will track the growth of a single instability wave
    % with x
    data=os_spatial_solver_x_sweep(re,om,nmodes,nxp,nyp,bl_prof);
    dimms=size(data);
 
    % Store results (spatial growth, i.e. imaginary alpha)
    for j=1:dimms(1)
        xloc(i,j)=data(j,1);
        growth(i,j)=data(j,2);
        cph(i,j)=data(j,3);
    end
    istart(i)=data(1,4);
    for j=1:dimms(1)
        nfact(i,j)=0;
    end

end

% Now loop over omega and x and compute n-factor
for i=1:nmax
    for j=istart(i):dimms(1)-1
        dx=(xloc(i,j+1)-xloc(i,j-1))/2;
        if cph(i,j)==0
            dn=0;
        else
            dn=dx*growth(i,j);
        end

        nfact(i,j)=nfact(i,j-1)+dn;
    end

end

% Plot data
for i=1:nmax

    clear lx ln lg
    for j=1:dimms(1)-1
        lx(j)=xloc(i,j);
    end
    %     for j=istart(i):dimms(1)-1
    for j=1:dimms(1)-1
        ln(j)=nfact(i,j);
        lg(j)=growth(i,j);
        lcph(j)=cph(i,j);
    end

    astr(i)=cellstr(sprintf('%d',omstore(i)));
    lx
    ln
    size(lx)
    size(ln)

    figure(1)
    subplot(1,1,1)
    plot(lx,lg)
    xlabel('x/c','FontSize',16)
    ylabel('imaginary wavenumber','FontSize',16)
    hold all
    figure(2)
    subplot(1,1,1)
    h(i)=plot(lx,ln);
    xlabel('x/c','FontSize',16)
    ylabel('ln(A/A0)','FontSize',16)
    hold all

    % figure(2)
    %
    % for i=1:(n+1)
    %     clear lcph
    %     for j=1:dimms(1)-1
    %         lcph(j)=cph(i,j);
    %     end
    %
    %     astr(i)=cellstr(sprintf('%d',omstore(i)));
    %     subplot(1,1,1)
    %
    %  h(i)=plot(lx,lcph);
    %     xlabel('x/c','FontSize',16)
    %     ylabel('cph','FontSize',16)
    %     hold all
    %     hold all
    %
    % end

    % save('PHASE_SPEED.dat', 'lcph', '-ASCII')
    % save('XC.dat', 'lx', '-ASCII')

    lnx=size(lx);

    %WRITE N FACTOR
    for k=1:lnx(2);
        datfile(1,k)=lx(k);
        datfile(2,k)=ln(k);
    end
    datfile;
    f=floor(omstore(i));
    d=round((omstore(i)-f)*100);

    f;
    d;

    title=sprintf('n_factor_omega_%d_%d.dat',f,d);
    title;

    fid = fopen(title,'w');

    % fprintf(fid,'%6.2f %12.8f\n',datfile(i,:));
    fprintf(fid,'%12.8f %12.8f\n',datfile);
    fclose(fid);


    %WRITE IMAG ALPHA
    for k=1:lnx(2);
        datfile(1,k)=lx(k);
        datfile(2,k)=lg(k);
    end
    datfile;
    f=floor(omstore(i));
    d=round((omstore(i)-f)*100);

    f;
    d;

    title=sprintf('a_i_omega_%d_%d.dat',f,d);
    title;

    fid = fopen(title,'w');

    % fprintf(fid,'%6.2f %12.8f\n',datfile(i,:));
    fprintf(fid,'%12.8f %12.8f\n',datfile);
    fclose(fid);


    %WRITE PHASE SPEED
    for k=1:lnx(2);
        datfile(1,k)=lx(k);
        datfile(2,k)=lcph(k);
    end
    datfile;
    f=floor(omstore(i));
    d=round((omstore(i)-f)*100);

    f;
    d;

    title=sprintf('c_ph_omega_%d_%d.dat',f,d);
    title;

    fid = fopen(title,'w');


    % fprintf(fid,'%6.2f %12.8f\n',datfile(i,:));
    fprintf(fid,'%12.8f %12.8f\n',datfile);
    fclose(fid);

end

max_n=0;
    omega_max=0;
    legend(h,astr,'Location','NW');
    for i=1:nmax;
        if nfact(i,dimms(1)-1) > max_n
            max_n=nfact(i,dimms(1)-1);
            omega_max=omstore(i);
        end
    end
    max_n
    omega_max

