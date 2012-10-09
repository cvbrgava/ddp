ord=10;
dia=ones(1,ord)*-2;
dia(1,ord)=dia(1,ord)+1;
A1=diag(dia,0);
A2=diag(ones(1,ord-1),1);
A3=diag(ones(1,ord-1),-1);
A=A1+A2+A3;
B=zeros(ord,1);
B(1,1)= B(1,1)+1;
C=zeros(1,ord);
C(1,ord)= C(1,ord)+1;
D=[0];
varVect=sym(zeros(1,ord));
for k=1:ord
    varVect(k)=sym(sprintf('v%d',k));
end

nonl=sym(zeros(ord,1));
nonl(1,1)= 2-exp(10*varVect(1))-exp(10*(varVect(1)-varVect(2)));
nonl(ord,1)=exp(10*(varVect(ord-1)-varVect(ord)))-1;
for k=2:ord-1
    nonl(k,1)= exp(10*(varVect(k-1)-varVect(k)))-exp(10*(varVect(k)-varVect(k+1)));
end


Anew=A*transpose(varVect)+nonl;
AJac=jacobian(Anew,varVect);
% AJac_subs=subs(AJac,varVect,zeros(1,ord));
% [b_act,a_act]=ss2tf(AJac_subs,B,C,D);
% H_act=tf(b_act,a_act);
% step(H_act)
% T=0:0.01:30;
% U=ones(size(T));
% [Yred,Xred]=lsim(AJac_subs,B,C,D,U,T);


%[Ynew1,Xnew1]=lsim(Anew,Bnew,Cnew,Dnew,U,T);
% plot(nonL10(:,1),nonL10(:,2),T,Y)%,T,Ynew,T,Ynew1)  
%figure(2)
%plot(T,Ynew) 

stTrj=[nonL1(:,2) nonL2(:,2) nonL3(:,2) nonL4(:,2) nonL5(:,2) nonL6(:,2) nonL7(:,2) nonL8(:,2) nonL9(:,2) nonL10(:,2)];
linPnts=[1 25 85 127 265]; linStates=zeros(ord,length(linPnts));
for index=1:length(linPnts)
    linStates(:,index)= transpose(stTrj(linPnts(index),:));
end
valinPnts=zeros(ord,length(linPnts));
Acurrent=zeros(ord,ord,length(linPnts));
for index=1:length(linPnts)
    Acurrent(:,:,index)=subs(AJac,varVect,transpose(linStates(:,index)));
    valinPnts(:,index)=subs(Anew,varVect,transpose(linStates(:,index)))-(Acurrent(:,:,index)*(linStates(:,index)));
    
end

red_ord=3;
X=ones(ord,red_ord,length(linPnts))*0;
X_new=ones(ord,red_ord,length(linPnts))*0;
%Arnoldi Method for basis evaluation
% 
for lin=1:length(linPnts)
    X(:,1,lin)=(inv(Acurrent(:,:,lin))*B)/norm(inv(Acurrent(:,:,lin))*B);
for index=2:red_ord
    basis=inv(Acurrent(:,:,lin))*X(:,index-1,lin);
    
    proj=zeros(ord,1);
    for x=1:index-1
        basis=basis-(X(:,x,lin)*transpose(X(:,x,lin))*basis);
    end
    
    unitBasis=basis/norm(basis);
    X(:,index,lin)=unitBasis;
end


%     X_new=ones(ord,red_ord,length(linPnts))*0;
      X_new(:,1,lin)=(inv(Acurrent(:,:,lin))*valinPnts(:,lin))/norm(inv(Acurrent(:,:,lin))*valinPnts(:,lin));
% %Arnoldi Method for basis evaluation
% % 
for index=2:red_ord
    basis=inv(Acurrent(:,:,lin))*X(:,index-1,lin);
    
    proj=zeros(ord,1);
    for x=1:index-1
        basis=basis-(X(:,x,lin)*transpose(X(:,x,lin))*basis);
    end
    
    unitBasis=basis/norm(basis);
    X_new(:,index,lin)=unitBasis;
end 
end
Xall=horzcat(X(:,:,1),X(:,:,2),X(:,:,3),X(:,:,4),X(:,:,5));
Xall_new=horzcat(X_new(:,:,2),X_new(:,:,3),X_new(:,:,4),X_new(:,:,5));
Complete=horzcat(Xall,Xall_new,linStates);
[eigenvec,eigenvals]=svd(Complete);
redBasis=[eigenvec(:,1) eigenvec(:,2) eigenvec(:,3)];
for index=1:length(linPnts)
    Acurrent_red(:,:,index)= transpose(redBasis)*Acurrent(:,:,index)*redBasis;
end
B_red=transpose(redBasis)*B;
linStates_new=(inv(eigenvec)*linStates);
linStates_red=linStates_new(1:3,:);
