%% If you use this algorithm, please cite
% Adeyemo, S., Bhattacharyya, D. Optimal nonlinear dynamic sparse model selection and...
% Bayesian parameter estimation for nonlinear systems. Comput. Chem. Eng. 180, 108502. https://doi.org/10.1016/j.compchemeng.2023.108502

%% Load data 
% load('data.mat');
% load('t.mat');
% load('u.mat');

% The rows of 'data' and 'u' are the measured output and input variables
% respectively while the columns are the time instants. 

%% Dimension of data
N = size(data,1);                                       % No. of output variablees
N2 = N^2;                       
[lengthu, lengthu2] = size(u);                          % Dimension of input data
lenreading = size(data,2)-1;

%% Range of input/output variables
umax=zeros(1,lengthu); umin=umax;
for i=1:lengthu
    umax(i) = max(u(i,:)) + 0.05*abs(max(u(i,:)));
    umin(i) = min(u(i,:)) - 0.05*abs(min(u(i,:)));
end

mxd=zeros(1,N); mnd=umax;
for i=1:size(data,1)
    mxd(i)=max(data(i,:)) + 0.05*abs(max(data(i,:)));
    mnd(i)=min(data(i,:)) - 0.05*abs(min(data(i,:)));
end


%% Normalization of input/output data
for i=1:size(data,1)
    data(i,:)=1 + 2*(data(i,:)-mnd(i))/(mxd(i)-mnd(i));
end
for i = 1:size(u,1)
    u(i,:) = 1+2*(u(i,:)-umin(i))/(umax(i)-umin(i));
end

%% Compute input transformations
[Utrans, iUY, i_UiUj] = InputTransformation(data,u,[],[],true);
Ntrans = size(Utrans,1)                                     % No. of transformations

UYindex1 = BasisInterpretation(lengthu,iUY,i_UiUj,[],[],[],false);


%% Selecting specific range(s) of data for ranking, model estimation and finetuning
rank_sel = 1:800; 
est_sel = 1:800;

data_rank = data(:,rank_sel);
u_rank = u(:,rank_sel);
utrans_rank = Utrans(:,rank_sel);

data_est = data(:,est_sel);
u_est = u(:,est_sel);
utrans_est = Utrans(:,est_sel);

%% Save selected data, inputs and their transformations
    % save('InputOutput','data1b','data2b','data3b','u1b','u2b','u3b','utrans1b','utrans2b','utrans3b')
    % clear data u Utrans data1 data2 data3 u1 u2 u3 utrans1 utrans2 utrans3 t1 t2 t3
    % clear data u Utrans

%% Ranking basis functions
save('data4rank','data_rank','u_rank','utrans_rank','UYindex1')
clear data_rank u_rank utrans_rank

gh=zeros(1,Ntrans,Ntrans);
for i=1:Ntrans
    gh(:,i,i)=1;
end

ind=zeros(1,Ntrans);
for i=1:Ntrans
    ind=gh(:,:,i); i
    [Ob(i,:), Am(i,:), Cm(i,:)] = RankBasis(ind,iUY,i_UiUj,UYindex1);
end
[Ob, idOb]=sort(Ob);
Utrans=Utrans(idOb,:);

%% 
utrans_est = utrans_est(idOb,:);
if Ntrans>85
    nut = 85;                                                   % No. of top rank Utrnas to be considered for BnB
    Ntrans = nut;
    utrans_est = utrans_est(1:nut,:);
else
    nut = Ntrans;
end

UYindex2 = [];
for i=1:numel(UYindex1)
    uyi = find(idOb==UYindex1(i));
    if uyi <= nut
        UYindex2 = [UYindex2 uyi];
    end
end

save('data4est','data_est','u_est','utrans_est','idOb','nut','UYindex2')
clear data_est u_est utrans_est u1 u2 u3 uu1 uu2 uu3 ut1 ut2 ut3

%% Branch and Bound for model selection
ndes = 10;                                                      % Number of desired hierarchically ranked models
Nworker = 8;                                                    % Number of workers
[AICcve, Amat, Cmat, Lmat, Obtrnd] = BranchNBound(N,ndes,Ntrans,iUY,i_UiUj,UYindex2,Nworker)

%% Save solution workspace
save('solution_workspace','N','AICcve','Amat','Cmat','est_sel','idOb','Lmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj','iUY','UYindex2')

%% Post processing of indices of transformed inputs for interpretation
% slnrank = 1;                                                    % Rank of model to be interpreted
% UYindex = BasisInterpretation(lengthu,iUY,i_UiUj,Obtrnd,idOb,slnrank,true);

%% Simulate model 1 and display training plots
slnrank = 1;
figure
plot(nonzeros(Obtrnd(slnrank,:)))
ylabel('Objective fcn (AICc+2Cve)')
xlabel('Number of Basis fcns')
title(['Trend of Objective Function for model ', num2str(slnrank)])

load('data4est.mat')
for i=1:N 
    data_est(i,:) = mnd(i) + 0.5*(data_est(i,:)-1)*(mxd(i)-mnd(i));
end
for i=1:lengthu 
    u_est(i,:) = umin(i) + 0.5*(u_est(i,:)-1)*(umax(i)-umin(i));
end
y = simulate_model_dynamic(u_est, data_est(:,1), slnrank);

for k=1:N
    figure
    plot(y(k,:),'Color',"#77AC30")
    hold on
    plot(data_est(k,:),'-.','Color',"#D95319")
    hold off
    legend('BIDSAM model','Data')
    title(['Training plot for variable y', num2str(k)])
    xlabel('time indices')
    ylabel(['y',num2str(k)])
end
