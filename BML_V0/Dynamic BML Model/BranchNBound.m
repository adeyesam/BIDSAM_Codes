function [Aic,Amat,Cmat,Lmat,Obtrnd]=BranchNBound(N,ndes,Ntrans,iUY,i_UiUj,UYindex,Nworker)
% ndes=7; mv=3; cv=6; Nworker=4;

% ndes is the number of solutions you desire, ie. 10 => returns the global top 10 solutions
% Nworker is the number of workers, ie processors, available to you to solve the problem

N2 = N^2; 
szlam = (N+N2)/2;
h1j = zeros(1,Ntrans);
timer = 0;
pct = 1e-5;

%h1j is a vector which you use to for your branching strategy
h1j=1:Ntrans;
flag1=true;
flag2=true;

%F is the fixed set of a node, C is the candidate set of a node, both vectors of logical 1s and 0s
F=false(1,Ntrans);
C=true(1,Ntrans);

%ASTB is your initial incumbent, ie perfection
ASTB=inf; AMATB=zeros(1,N2); CMATB=zeros(1,N*Ntrans); LMATB=zeros(1,szlam); OBTRNDB = zeros(1,Ntrans);


ita=0;
Aic=Inf(ndes,1);
Amat=zeros(ndes,N2);
Cmat=zeros(ndes,N*Ntrans);
Lmat=zeros(ndes,szlam);
% sset=zeros(ndes,mv);
Obtrnd=zeros(ndes,Ntrans);
dumbsolve();

    function dumbsolve()
        while ~isempty(F)
            ita=ita+1;
            Abound=max(Aic);
            %Pruning of supernodes
            if ita>2
                [bumb,~]=size(F);
                discard=false(bumb,1);
                L1=ASTB>Abound;
                discard=find(L1);
                if any(discard)
                    F(discard,:)=[];
                    C(discard,:)=[];
                    ASTB(discard,:)=[];
                    AMATB(discard,:)=[];
                    CMATB(discard,:)=[];
                    LMATB(discard,:)=[];
                    OBTRNDB(discard,:)=[];
                end
            end
            if ~isempty(F)    
                %Tracks Progress
                FC=[];
                CC=[];
                ASTBT=[];
                AMATBT=[];
                CMATBT=[];
                LMATBT=[];
                OBTRNDBT=[];
                [NumberofNodes,~]=size(F)
                t=0;
                if NumberofNodes>28000
                    flagconst=true;
                else
                    flagconst=false;
                end
            end   
            while NumberofNodes<Nworker && ~isempty(F)
                [~,idc]=max(ASTB);               
                f1=F(idc,:);
                c1=C(idc,:);
                f2=f1;
                h1k=h1j;
                h1k(~c1)=Inf;
                [~,id]=min(h1k);
                f2(id)=true;
                c1(id)=false;
                F(idc,:)=[];
                C(idc,:)=[];
                [a1,amat1,cmat1,lmat1]=feval(@Gc,f2,Ntrans,ASTB(idc,:),AMATB(idc,:),CMATB(idc,:),LMATB(idc,:),N,iUY,i_UiUj,UYindex,true,false);
                [a2,amat2,cmat2,lmat2]=feval(@Gc,f1,Ntrans,ASTB(idc,:),AMATB(idc,:),CMATB(idc,:),LMATB(idc,:),N,iUY,i_UiUj,UYindex,false,false);
                
                obtrndp = OBTRNDB(idc,:);
                percent1=(ASTB(idc,:)-a1)/abs(ASTB(idc,:));

                if (ita==1)||(sum(f2)==0) && a1 <= ASTB(idc,:)
                    ASTB = [ASTB; a1];
                    F=[F ; f2];
                    C=[C ; c1];
                    AMATB=[AMATB ; amat1];
                    CMATB=[CMATB ; cmat1];
                    LMATB=[LMATB ; lmat1];
                    obtrnd1=obtrndp; obtrnd1(:,id)=a1;
                    OBTRNDB=[OBTRNDB ; obtrnd1]; 

                elseif sum(c1)==0 || percent1<= pct
                    Aic=[Aic ; a1];
                    Amat=[Amat ; amat1];
                    Cmat=[Cmat ; cmat1];
                    Lmat=[Lmat ; lmat1];
                    obtrnd1=obtrndp; obtrnd1(:,id)=a1;
                    Obtrnd=[Obtrnd ; obtrnd1];

                elseif a1 < ASTB(idc,:) && percent1> pct
                    ASTB = [ASTB; a1];
                    F=[F ; f2];
                    C=[C ; c1];%CMATB
                    AMATB=[AMATB ; amat1];
                    CMATB=[CMATB ; cmat1];
                    LMATB=[LMATB ; lmat1];
                    obtrnd1=obtrndp; obtrnd1(:,id)=a1;
                    OBTRNDB=[OBTRNDB ; obtrnd1];                    
                end

                if sum(c1)>0
                    ASTB = [ASTB; a2];
                    F=[F ; f1];
                    C=[C ; c1];
                    AMATB=[AMATB ; amat2];
                    CMATB=[CMATB ; cmat2]; 
                    LMATB=[LMATB ; lmat2]; 
                    obtrnd2=obtrndp; %obtrnd(:,find(f1,1,'last'))=a2;
                    OBTRNDB=[OBTRNDB ; obtrnd2];
                 end
                
                if isempty(F) && sum(c1)==0 %&& ASTB(idc,:)>max(Aic)
                    Aic=[Aic ; a2];
                    Amat=[Amat ; amat2];
                    Cmat=[Cmat ; cmat2];
                    Lmat=[Lmat ; lmat2];
                    Obtrnd=[Obtrnd ; obtrnd2];
                end
                
                ASTB(idc,:)=[];
                AMATB(idc,:)=[];
                CMATB(idc,:)=[];
                LMATB(idc,:)=[];
                OBTRNDB(idc,:)=[];
                
                [NumberofNodes,~]=size(F);
            end
            if ~isempty(F)  
                iterat=inf;
                if ~flagconst
                    if NumberofNodes<100*Nworker
                        FC=F;
                        F=[];
                        CC=C;
                        C=[];
                        ASTBT=ASTB;
                        ASTB=[];
                        AMATBT=AMATB;
                        AMATB=[];
                        CMATBT=CMATB;
                        CMATB=[];
                        LMATBT=LMATB;
                        LMATB=[];
                        OBTRNDBT=OBTRNDB;
                        OBTRNDB=[];
                    else
                        shf=randperm(NumberofNodes,100*Nworker);
                        FC=F(shf,:);
                        F(shf,:)=[];
                        CC=C(shf,:);
                        C(shf,:)=[];
                        ASTBT=ASTB(shf,:);
                        ASTB(shf,:)=[];
                        AMATBT=AMATB(shf,:);
                        AMATB(shf,:)=[];
                        CMATBT=CMATB(shf,:);
                        CMATB(shf,:)=[];
                        LMATBT=LMATB(shf,:);
                        LMATB(shf,:)=[];
                        OBTRNDBT=OBTRNDB(shf,:);
                        OBTRNDB(shf,:)=[];
                    end
                else
                    [~,shf]=sort(sum(C,2));
                    F=F(shf,:);
                    C=C(shf,:);
                    ASTB=ASTB(shf,:);
                    AMATB=AMATB(shf,:);
                    CMATB=CMATB(shf,:);
                    LMATB=LMATB(shf,:);
                    OBTRNDB=OBTRNDB(shf,:);
                    FC=F(1:100*Nworker,:);
                    F(1:100*Nworker,:)=[];
                    CC=C(1:100*Nworker,:);
                    C(1:100*Nworker,:)=[];
                    ASTBT=ASTB(1:100*Nworker,:);
                    ASTB(1:100*Nworker,:)=[];
                    AMATBT=AMATB(1:100*Nworker,:);
                    AMATB(1:100*Nworker,:)=[];
                    CMATBT=CMATB(1:100*Nworker,:);
                    CMATB(1:100*Nworker,:)=[];
                    LMATBT=LMATB(1:100*Nworker,:);
                    LMATB(1:100*Nworker,:)=[];
                    OBTRNDBT=OBTRNDB(1:100*Nworker,:);
                    OBTRNDB(1:100*Nworker,:)=[];
                end
                [NodestoWorker,~]=size(FC)
                Li=randperm(NodestoWorker); 
                FC=FC(Li,:);
                CC=CC(Li,:);  
                ASTBT=ASTBT(Li,:);
                AMATBT=AMATBT(Li,:);
                CMATBT=CMATBT(Li,:);
                LMATBT=LMATBT(Li,:);
                OBTRNDBT=OBTRNDBT(Li,:);
                divid=floor(NodestoWorker/Nworker);
                for tt=1:Nworker
                    if tt==1
                        t=1;
                    else
                        t=t+divid;
                    end
                    if tt==Nworker
                        problems(tt).F=FC(t:NodestoWorker,:);
                        problems(tt).C=CC(t:NodestoWorker,:);
                        problems(tt).A=ASTBT(t:NodestoWorker,:);
                        problems(tt).AMATB=AMATBT(t:NodestoWorker,:);
                        problems(tt).CMATB=CMATBT(t:NodestoWorker,:);
                        problems(tt).LMATB=LMATBT(t:NodestoWorker,:);
                        problems(tt).OBTRNDB=OBTRNDBT(t:NodestoWorker,:);
                    else
                        problems(tt).F=FC(t:(divid*tt),:);
                        problems(tt).C=CC(t:(divid*tt),:);
                        problems(tt).A=ASTBT(t:(divid*tt),:);
                        problems(tt).AMATB=AMATBT(t:(divid*tt),:);
                        problems(tt).CMATB=CMATBT(t:(divid*tt),:);
                        problems(tt).LMATB=LMATBT(t:(divid*tt),:);
                        problems(tt).OBTRNDB=OBTRNDBT(t:(divid*tt),:);
                    end
                end

                %This section is solved in parallel
                parfor z=1:Nworker
                  if ~isempty(problems(z).F)  
                    FT=problems(z).F
                    CT=problems(z).C;
                    AT=Inf;
                    AST=problems(z).A;
                    AMATST=problems(z).AMATB;
                    CMATST=problems(z).CMATB;
                    LMATST=problems(z).LMATB;
                    OBTRNDST=problems(z).OBTRNDB;
                    AMATT=zeros(1,N2);
                    CMATT=zeros(1,N*Ntrans);
                    LMATT=zeros(1,szlam);
                    OBTRNDT=zeros(1,Ntrans);
                    DUMBA=Aic;
                    timer=0;
                    termin=true;
                    toccers=0;
                    while termin
                        tic
                        %chooses which node to branch
                        timer=timer+1
                        if flagconst
                            [~,idn]=min(sum(CT,2));
                        else
                            if toccers<90
                                [~,idn]=min(AST);
                            else
                                [~,idn]=min(sum(CT,2));
                            end
                        end
                        f1=FT(idn,:)
                        c1=CT(idn,:)
                        AS=AST(idn,:);
                        AMATS=AMATST(idn,:);
                        CMATS=CMATST(idn,:);
                        LMATS=LMATST(idn,:);
                        OBTRNDS=OBTRNDST(idn,:);
                        f2=f1;

                        %chooses which variable to branch
                        h1k=h1j;
                        h1k(~c1)=Inf
                        [~,id]=min(h1k);
                        f2(id)=true
                        c1(id)=false
                        c2=c1;
                        OBTRNDST_idn=OBTRNDST(idn,:);

                        %removing node to be branched
                        FT(idn,:)=[];
                        CT(idn,:)=[];
                        AST(idn,:)=[];
                        AMATST(idn,:)=[];
                        CMATST(idn,:)=[];
                        LMATST(idn,:)=[];
                        OBTRNDST(idn,:)=[];
                        flag1=true;
                        flag2=true;

                        %checking whether node1 is terminal
                        cc=sum(c1);
                        if flag1 && ~isempty(OBTRNDST_idn)
                          [a,amat,cmat,lmat]=feval(@Gc,f2,Ntrans,0,zeros(1,N2),zeros(1,N*Ntrans),zeros(1,szlam),N,iUY,i_UiUj,UYindex,true,true);
                          Abound=min([max(AT) max(DUMBA)]);
                          %if a<AS %Abound
                                percent=(AS-a)/abs(AS)
                            if (cc==0 || percent<= pct) &&(a<AS)
                                AT=[a ; AT];
                                AMATT=[amat ; AMATT];
                                CMATT=[cmat ; CMATT];
                                LMATT=[lmat ; LMATT];
                                DUMBA=[a ; DUMBA];
                                obtrnd=OBTRNDS; obtrnd(:,id)=a;
                                OBTRNDT=[obtrnd ; OBTRNDT];

                                if size(AT,1)>ndes
                                    [~,IDDD]=max(AT);
                                    AT(IDDD,:)=[];
                                    AMATT(IDDD,:)=[];
                                    CMATT(IDDD,:)=[];
                                    LMATT(IDDD,:)=[];
    %                                 ssetT(IDDD,:)=[];
                                    OBTRNDT(IDDD,:)=[];
                                end

                                flag1=false;                           
                            end
                        end

                        %checking whether node2 is terminal
                        cc=sum(c2)
                       
                        if flag2 && cc==0
%                                 if cc==0
                                AT=[AS ; AT];
                                AMATT=[AMATS ; AMATT];
                                CMATT=[CMATS ; CMATT];
                                LMATT=[LMATS ; LMATT];
                                DUMBA=[AS ; DUMBA];
                                OBTRNDT=[OBTRNDS ; OBTRNDT];
                                
                                if size(AT,1)>ndes
                                    [~,IDDD]=max(AT);
                                    AT(IDDD,:)=[];
                                    AMATT(IDDD,:)=[];
                                    CMATT(IDDD,:)=[];
                                    LMATT(IDDD,:)=[];
                                    OBTRNDT(IDDD,:)=[];
                                end
                                
                                flag2=false;
                        end
                        
                        %Pruning of Branch 1
                        if (flag1 && cc>0)
                            [a,amat,cmat,lmat]=feval(@Gc,f2,Ntrans,AS,AMATS,CMATS,LMATS,N,iUY,i_UiUj,UYindex,true,false);
                            Abound=min([max(AT) max(DUMBA)]);
                            percent=(AS-a)/abs(AS)
                            if a<AS && percent>pct 
                                AST=[AST ; a];
                                FT=[FT ; f2];
                                CT=[CT ; c1];
                                AMATST=[AMATST ; amat];
                                CMATST=[CMATST ; cmat];
                                LMATST=[LMATST ; lmat];

                                obtrnd=OBTRNDS; obtrnd(:,id)=a;
                                OBTRNDST=[OBTRNDST ; obtrnd];
                            else
                                flag1=false;
                            end
                        end

                        %Pruning of Branch 2
                        if (flag2 && cc>0)
                                AST=[AST ; AS];                            
                                FT=[FT ; f1];
                                CT=[CT ; c2];
                                AMATST=[AMATST ; AMATS];
                                CMATST=[CMATST ; CMATS];
                                LMATST=[LMATST ; LMATS];

                                obtrnd=OBTRNDS; 
                                OBTRNDST=[OBTRNDST ; obtrnd];
                        else
                                flag2=false;
                        end
                        dtt=toc;
                        toccers=toccers+dtt;
                        if timer>iterat || isempty(FT) || toccers>300
                            termin=false;
                        end
                   end
                    %Storing existing nodes and solutions from iteration
                    [ttte,~]=size(FT);
                    if ttte>0
                        F=[F ; FT];
                        C=[C ; CT];
                        ASTB=[ASTB ; AST];
                        AMATB=[AMATB ; AMATST];
                        CMATB=[CMATB ; CMATST];
                        LMATB=[LMATB ; LMATST];
                        OBTRNDB=[OBTRNDB ; OBTRNDST];
                    end
                    SOLUTIONS(z).D2=AT;
                    SOLUTIONS(z).D3=AMATT;
                    SOLUTIONS(z).D4=CMATT;
                    SOLUTIONS(z).D5=LMATT;
                    SOLUTIONS(z).D6=OBTRNDT;
                  end  
                end

                %Post-processing of solutions
                [~,st]=size(SOLUTIONS)
                for i=1:st
                    Aic=[Aic ; SOLUTIONS(i).D2];
                    Amat=[Amat ; SOLUTIONS(i).D3];
                    Cmat=[Cmat ; SOLUTIONS(i).D4];
                    Lmat=[Lmat ; SOLUTIONS(i).D5];
                    Obtrnd=[Obtrnd ; SOLUTIONS(i).D6];
                end
            end
            if size(Aic,1)>ndes
                [Aic,idpd]=sort(Aic);
                Amat=Amat(idpd,:);
                Cmat=Cmat(idpd,:);
                Lmat=Lmat(idpd,:);
                Obtrnd=Obtrnd(idpd,:);
                idsss=size(Aic);
                Aic((ndes+1):idsss,:)=[];
                Amat((ndes+1):idsss,:)=[];
                Cmat((ndes+1):idsss,:)=[];
                Lmat((ndes+1):idsss,:)=[];
                Obtrnd((ndes+1):idsss,:)=[];
            end
            SOLUTIONS=[];
        end
    end


end

%This is your objective function
function [a,amat,cmat,lmat]=Gc(f,cv,as,amats,cmats,lmats,n,iUY,i_UiUj,UYindex,flag,flag2)
szlam=n^2 + (n-n^2)/2;
k=cv-sum(f);
warning('off', 'all')
if flag2
        [a, amat, cmat, lmat]=ModelEstimation(f,iUY,i_UiUj,UYindex,amats,cmats,lmats);
else
    if k==cv
        a=inf; amat=zeros(1,n*n); cmat=zeros(1,n*cv); lmat=zeros(1,szlam);
    elseif flag
        [a, amat, cmat, lmat]=ModelEstimation(f,iUY,i_UiUj,UYindex,amats,cmats,lmats);
    else
        a=as; amat=amats; cmat=cmats; lmat=lmats;
    end
end
end