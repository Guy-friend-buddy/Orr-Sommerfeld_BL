function bl_dat=bl_reader

fid = fopen('C:\Users\pj2v20\Downloads\Orr-Sommerfeld\Orr-Sommerfeld\Convective Instability Growth\bl_data.dat','r');
source = fscanf(fid,'%g');
status = fclose(fid);

nxp=source(1);
nyp=source(2);
size(source);
for i=1:nxp;
    xloc(i)=source(i+2);
end

for i=1:nxp;
    
    for j=1:nyp;
        
        u(i,j)=source(2+nxp+nyp*(i-1)+j);
        y(i,j)=source(2+nxp+nyp*(i-1+nxp)+j);
    end

end


    
    for i=1:nxp
        for j=1:nyp/3
            u2(i,j)=u(i,j);
            y2(i,j)=y(i,j);
        end
    end
    


for i=1:10:nxp
    
        subplot(1,2,1)
%         figure(1)       
        plot(u2(i,:),y2(i,:))
hold on
end

for i=1:10:nxp
    
        subplot(1,2,2)
%        figure(2)
        plot(u(i,:),y(i,:))

end

hold off

for j=1:nyp
    for i=1:nxp
    bl_dat(i,j,1)=u(i,j);
    bl_dat(i,j,2)=y(i,j);
    end
end

j=nyp+1;
for i=1:nxp
bl_dat(i,j,1)=xloc(i);
end


    
    
    