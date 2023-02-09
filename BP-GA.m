clc
clear
file_name = dir(fullfile('E:\secdata\','*20210911sec.sec'));    
File_name = dir(fullfile('E:\secdata\','*20210912sec.sec'));    
LengthFiles = length(file_name);
format long 
for p=1:LengthFiles    
  fid=fopen(strcat('E:\secdata\',file_name(p).name),'r');     
  Fid=fopen(strcat('E:\secdata\',File_name(p).name),'r');     
  for ii=1:13
    line1=fgetl(fid)
    line2=fgetl(Fid)
  end
  k=1;
  while ~feof(fid)
      line1=fgetl(fid)
      line2=fgetl(Fid)
      D(k,p)=str2double(line1(33:40));
      H(k,p)=str2double(line1(43:50));
      Z(k,p)=str2double(line1(53:60));
      DD(k,p)=str2double(line2(33:40));
      HH(k,p)=str2double(line2(43:50));
      ZZ(k,p)=str2double(line2(53:60));      
      k=k+1;
  end
end
fclose(Fid);
fclose(fid);
D(:,[1,2])=D(:,[2,1]);
H(:,[1,2])=H(:,[2,1]);
DD(:,[1,2])=DD(:,[2,1]);
Z(:,[1,2])=Z(:,[2,1]);
HH(:,[1,2])=HH(:,[2,1]);
ZZ(:,[1,2])=ZZ(:,[2,1]);


P=D(:,2:5);    
P=P';
T=D(:,1);       
T=T';
P1=DD(:,2:5);      
P1=P1';
T1=DD(:,1);         
T1=T1';


inputnum=4;
hiddennum=5;
outputnum=1;
inputnum=4;
hiddennum=5;
outputnum=1;


[pn,pp]=mapminmax(P,0,1);
[tn,tt]=mapminmax(T,0,1);
[pn1,pp1]=mapminmax(P1,0,1);
[tn1,tt1]=mapminmax(T1,0,1);
net=newff(minmax(pn),[5,1],{'logsig','purelin','traingdm'});    


maxgen=15;                         
sizepop=15;                        
pcross=[0.9];                     
pmutation=[0.01];                   

numsum=inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum;

lenchrom=ones(1,numsum);        
bound=[-3*ones(numsum,1) 3*ones(numsum,1)];   

%--------------------------------------------------------------------------------------------------------------
individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]);  
avgfitness=[];                    
bestfitness=[];                     
bestchrom=[];                      

for i=1:sizepop
   
    individuals.chrom(i,:)=Code(lenchrom,bound);    
    x=individuals.chrom(i,:);
   
    individuals.fitness(i)=fun(x,inputnum,hiddennum,outputnum,net,pn,tn);   %染色体的适应度
end
FitRecord=[];

[bestfitness bestindex]=min(individuals.fitness);
bestchrom=individuals.chrom(bestindex,:); 
avgfitness=sum(individuals.fitness)/sizepop; 

trace=[avgfitness bestfitness]; 


for i=1:maxgen
    i
 
    individuals=Select(individuals,sizepop); 
    avgfitness=sum(individuals.fitness)/sizepop;

    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);

    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,i,maxgen,bound);
    
 
    for j=1:sizepop
        x=individuals.chrom(j,:); 
        individuals.fitness(j)=fun(x,inputnum,hiddennum,outputnum,net,pn,tn);   
    end
    
 
    [newbestfitness,newbestindex]=min(individuals.fitness);
    [worestfitness,worestindex]=max(individuals.fitness);

    if bestfitness>newbestfitness
        bestfitness=newbestfitness;
        bestchrom=individuals.chrom(newbestindex,:);
    end
    individuals.chrom(worestindex,:)=bestchrom;
    individuals.fitness(worestindex)=bestfitness;
    
    avgfitness=sum(individuals.fitness)/sizepop;
    
    trace=[trace;avgfitness bestfitness]; 
    FitRecord=[FitRecord;individuals.fitness];
end


figure(1)
[r c]=size(trace);
plot([1:r]',trace(:,2),'b--');
title(['fitness' 'No＝' num2str(maxgen)]);
xlabel('No');ylabel('fitness');
legend('average','best');
disp('finess');



w1=x(1:inputnum*hiddennum);
B1=x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);
w2=x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);
B2=x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum);

net.iw{1,1}=reshape(w1,hiddennum,inputnum);
net.lw{2,1}=reshape(w2,outputnum,hiddennum);
net.b{1}=reshape(B1,hiddennum,1);
net.b{2}=B2;

net.trainParam.show=50;                                                         

net.trainParam.Ir=0.05;                                                                  
net.trainParam.mc=0.9;
net.trainParam.epochs=1000;                                                        
net.trainParam.goal=0.001;                                                          
[net,tr]=train(net,pn,tn);
A=sim(net,pn);       
Q = mapminmax('reverse',A,tt)   
V=nanmean(abs(T-Q));  
save quanqueHnet
A1=sim(net,pn1);
Q1 = mapminmax('reverse',A1,tt1)   

V1=nanmean(abs(T1-Q1));      

Q1 =Q1';

A=reshape(T1,1800,48);
B=reshape(Q1,1800,48);
parfor i=1:48
   det1(i)=A(1,i)-B(1,i);
   det2(i)=A(1800,i)-B(1800,i);
   det(i)=(det2(i)-det1(i))/1800;
      for j=1:1800
         rebnew(j,i)=j*det(i)+det1(i)+B(j,i)
   end
end
C=reshape(rebnew,86400,1);
Error=nanmean(abs(T1'-C));
fprintf('\nthe BP error is%9.4f\n',V);
fprintf('\nthe BP error is%9.4f\n',V1);
fprintf('\nthe rebuilt error is%9.4f\n',Error);
figure(2)
plot(C,'r');
hold on
plot(T1,'blue')



