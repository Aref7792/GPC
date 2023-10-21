# GPC
#Coded GPC algorithm presented in chapter 6 of "Model Predictive Control in the Process Industry" by By Eduardo F. Camacho, Carlos A. Bordons

clc;
clear;
close all;


%% system information

A={[1 -1.8586 .8636];[1 -2.7634 2.5452 -.7813];[1 -2.6091 2.2663 -.6553]};

B={[0 0 0 0 0 0 .0802 .1562 -.2163] ,[0 0 0 0 0 0 0 .1142 -.1054],[0 0 0 0 0 0 .1164 .2267 -.3140];[0 0 0 0 .2113 -.1858 -.1949 .1719],[0 0 0 .1875 -.1613 -.1750 .1515],[0 0 0 .1704 .1696 -.7567 .4199];[0 0 0 0 0 .5 -.8617 .3699 ],[0 0 0 0 0  0.1964    -.1454   -.1774    .1347],[1.367 -2.4591 1.1057]};

n=3; %num. of outputs

m=3; %num. of inputs

Np=30;
Nu=5;

%% calculating Ahat

for i=1:n
    Ahat{i,1}=conv([1 -1],A{i,1});
end

%% calculating E&F

for ii=1:n
    E{1,ii}=1;
    for k=2:Np+1
        for i=2:k
            e(1)=1;
            for j=1:i-1
                if j+1>numel(Ahat{ii,1})
                    aa=0;
                else
                    aa=Ahat{ii,1}(j+1);
                end
                s(j)=-aa*e(i-j);
            end
            e(i)=sum(s);
            clear s;
        end
        E{k,ii}=e;
        clear e
    end
end



for ii=1:n
    siz=size(A{ii,1});
    siz=siz(2);
    for k=1:Np+1
        for j=1:siz
            for i=1:k+j-1
                if i+1>numel(Ahat{ii,1})
                    aaa=0;
                else
                    aaa=Ahat{ii,1}(i+1);
                end
                if k+j-i>numel(E{k,ii})
                    ee=0;
                else
                    ee=E{k,ii}(k+j-i);
                end
                ss(i)=-aaa*ee;
            end
            f(j)=sum(ss);
            clear ss
        end
        F{k,ii}=f;
        clear f
    end
end

%% calculating input delays

for i=1:n
    for j=1:m
        d=1;
        
        num=numel(B{i,j});
        
        for ii=1:num
            if B{i,j}(ii)==0
                d=d+1;
            else
                break ;
            end
        end
        d=d-1;
        D(i,j)=d;
        clear d
    end
end

BB=B;
for i=1:n
    del=min(D(i,:));
    for j=1:m
        BB{i,j}(1:del)=[];
    end
end

%% calculating GG

for i=1:n
    for j=1:m
        for k=1:Np+1
            GG{i,1}{k,j}=conv(BB{i,j},E{k,i});
        end
    end
end

%% calculating G&Gp


for i=1:n
    N1(i)=min(D(i,:))+1;
    for j=1:m
        for k=N1(i):Np
            G{i,1}{k-N1(i)+1,j}=fliplr(GG{i,1}{k,j}(1:k-N1(i)+1));
            Gp{i,1}{k-N1(i)+1,j}=GG{i,1}{k,j}(k-N1(i)+2:end);
            
        end
    end
end

for i=1:n
    N1(i)=min(D(i,:))+1;
    for k=N1(i):Np
        FF{i,1}{k-N1(i)+1,1}=F{k,i};
    end
end
        

for i=1:n
    nr=size(G{i,1});
    nr=nr(1);
    for j=1:m
        ne(i)=numel(G{i,1}{nr,j});
    end
end


for i=1:n
    nr=size(G{i,1});
    nr=nr(1);
    for j=1:m
        ne=max(ne);
        for k=1:nr
            for ii=1:ne
                if ii>numel(G{i,1}{k,j})
                    GP{i,1}{k,j}(ii)=0;
                else
                    GP{i,1}{k,j}(ii)=G{i,1}{k,j}(ii);
                end
            end
        end
    end
end

for i=1:n
    nr=size(GP{i,1});
    nr=nr(1);
    for j=1:m
        for k=1:nr
            GN{i,1}{k,j}=GP{i,1}{k,j}(1:Nu);
        end
    end
end

for i=1:n
    GNN{i,1}=cell2mat(GN{i,1});
end

GNN=cell2mat(GNN);

 %% main loop
 
 nq=size(GNN);
 
 Q=.3*eye(nq(2));
 
 R=1*eye(nq(1));
 
 gain=(GNN'*R*GNN+Q)^-1*GNN'*R;
 
 for i=1:n
     for j=1:m
         neb(i,j)=numel(B{i,j});
     end
 end
 
 for i=1:n
     nea(i,1)=numel(A{i,1});
 end
 
 nebm=max(max(neb));
 
 neam=max(nea);
 
 nin=max(nebm,neam);
 
 y=zeros(nin,n);
 
 u=zeros(nin,m);
 
 du=zeros(nin,m);
 
 t=1:nin;
 
 for i=nin+1:100
     
     if i<51
         sp=[.5 .3 .1];
     else
         sp=[.4 .3 .1];
     end
     
    for j=1:n
        for k=1:m
            TB(k)=B{j,k}*u(i-1:-1:i-neb(j,k),k);
        end
        STB=sum(TB);
        y(i,j)=-A{j,1}(2:end)*y(i-1:-1:i-nea(j)+1,j)+STB;
    end
    
    for ii=1:n
        nk=size(Gp{ii,1});
        nk=nk(1);
        for j=1:m
            for k=1:nk
                nf=numel(Gp{ii,1}{k,j});
                ff{ii,1}{k,j}=Gp{ii,1}{k,j}*du(i-1:-1:i-nf,j);
            end
        end
    end
    
    for ii=1:n
        nk=size(Gp{ii,1});
        nk=nk(1);
        for k=1:nk
            for j=1:m
                sf(j)=ff{ii,1}{k,j};
            end
            fff{ii,1}{k,1}=sum(sf);
        end
    end
    
    for ii=1:n
        nk=size(FF{ii,1});
        nk=nk(1);
        for k=1:nk
            nF=numel(FF{ii,1}{k,1});
            fy{ii,1}{k,1}=FF{ii,1}{k,1}*y(i:-1:i-nF+1,ii);
        end
    end
    
    for ii=1:n
        fff{ii,1}=cell2mat(fff{ii,1});
        fy{ii,1}=cell2mat(fy{ii,1});
    end
    

    FFF=cell2mat(fff);
    FY=cell2mat(fy);
    
    clear fff
    clear fy
    
    Ff=FFF+FY;

    
    for ii=1:n
        nk=size(Gp{ii,1});
        nk=nk(1);
        W{ii,1}=repmat(sp(ii),nk,1);
    end
    
    WW=cell2mat(W);
    
    
    
    for ii=1:m
        nnu=(ii-1)*Nu+1;
        du(i,ii)=gain(nnu,:)*(WW-Ff);
    end
    
   
   
    
    for ii=1:m
        u(i,ii)=u(i-1,ii)+du(i,ii);
    end
    

    
    t(i)=i;
 end
 

plot(t,y);
grid on
%legend('u1','u2','u3');
