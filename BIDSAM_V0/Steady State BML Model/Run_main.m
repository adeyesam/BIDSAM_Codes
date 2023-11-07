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
lenreadpone = lenreading+1;

%% Range of input/output variables
umax=zeros(1,lengthu); umin=umax;
for i=1:lengthu
    umax(i) = max(u(i,:)) + 0.15*abs(max(u(i,:)));
    umin(i) = min(u(i,:)) - 0.15*abs(min(u(i,:)));
end
for i=1:size(data,1)
    mxd(i)=max(data(i,:)) + 0.15*abs(max(data(i,:)));
    mnd(i)=min(data(i,:)) - 0.15*abs(min(data(i,:)));
end


%% Normalization of input/output data
for i=1:size(data,1)
    data(i,:)=1 + 2*(data(i,:)-mnd(i))/(mxd(i)-mnd(i));
end
for i = 1:size(u,1)
    u(i,:) = 1+2*(u(i,:)-umin(i))/(umax(i)-umin(i));
end

%% Compute input transformations
[Utrans, i_UiUj] = InputTransformation(data,u,[],true);
Ntrans = size(Utrans,1)                                     % No. of transformations
BasisInterpretation(lengthu,i_UiUj,[],[],[],false);

%% Selecting specific data points for ranking and estimation
rank_sel = randsample(lenreadpone,floor(0.3*lenreadpone));
est_sel = randsample(lenreadpone,floor(0.7*lenreadpone));

data_rank=data(:,rank_sel);
u_rank=u(:,rank_sel);
utrans_rank=Utrans(:,rank_sel);

data_est=data(:,est_sel);
u_est=u(:,est_sel);
utrans_est=Utrans(:,est_sel);

%% Ranking basis functions
save('data4rank','data_rank','u_rank','utrans_rank')
clear data_rank u_rank utrans_rank

gh=zeros(1,Ntrans,Ntrans);
for i=1:Ntrans
    gh(:,i,i)=1;
end

ind=zeros(1,Ntrans);
for i=1:Ntrans
    ind=gh(:,:,i); i
    [Ob(i,:), Cm(i,:)] = RankBasis(ind);
end
[Ob, idOb]=sort(Ob);
Utrans=Utrans(idOb,:);

%%
utrans_est = utrans_est(idOb,:);
if Ntrans>85
    nut = 85;                                               % No. of top rank Utrnas to be considered for BnB
    Ntrans = nut;
    utrans_est = utrans_est(1:nut,:);
else
    nut = Ntrans;
end

save('data4est','data_est','u_est','utrans_est','idOb','nut')
clear data_est u_est utrans_est 

%% Branch and Bound for model selection
ndes = 8;                                                   % Number of desired hierarchically ranked models
Nworker = 8;                                                % Number of workers
[AICcve, Cmat, Lmat, Obtrnd] = BranchNBound(N,ndes,Ntrans,Nworker)

figure
plot(nonzeros(Obtrnd(1,:)))
ylabel('Objective fcn (AICc)')
xlabel('Number of Basis fcns')

%% Save solution workspace
save('solution_workspace','N','AICcve','Cmat','est_sel','idOb','Lmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj')

%% Post processing of indices of transformed inputs for interpretation
% slnrank = 1;                                                % Rank of model to be interpreted
% BasisInterpretation(lengthu,i_UiUj,Obtrnd,idOb,slnrank,true);

%% Simulate model 2 and display training plots
slnrank = 2;
load('data4est.mat')
for i=1:N 
    data_est(i,:) = mnd(i) + 0.5*(data_est(i,:)-1)*(mxd(i)-mnd(i));
end
for i=1:lengthu 
    u_est(i,:) = umin(i) + 0.5*(u_est(i,:)-1)*(umax(i)-umin(i));
end
y = simulate_model_steady(u_est, slnrank);

for k=1:N
    figure
    plot(data_est(k,:),y(k,:),'+','Color',"#77AC30")
    hold on
    plot(data_est(k,:),data_est(k,:),'-','Color',"#D95319")
    hold off
    % legend('BIDSAM model','Data')
    title(['Training plot for variable y', num2str(k)])
    xlabel('Data')
    ylabel('BIDSAM Model')
end
