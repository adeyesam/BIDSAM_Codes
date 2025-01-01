%% The file loads and simulates hierarchically ranked models (based on AICc+2Cve values)
% from the branch and bound algorithm. 'solnrank' is the rank of model we desire to simulate.
% The rows of 'inputs' are the input variables while the columns are the time instants


function y = simulate_model_dynamic(inputs, y_initial, solnrank)

load('solution_workspace')
u = inputs;
[lengthu, lengthu2] = size(u);
N2 = N^2;

% Scaling inputs and initial outputs
for i = 1:lengthu
    u(i,:) = 1+2*(u(i,:)-umin(i))/(umax(i)-umin(i));
end
for i=1:N
    y_initial(i,:)=1 + 2*(y_initial(i,:)-mnd(i))/(mxd(i)-mnd(i));
end

XXQ = y_initial;
XQ(:,1) = XXQ;


for jj=1:lengthu2
    
    % Compute input transformations
    if lengthu>1
        i_UiUj = nchoosek(1:lengthu,2);                     % Indices for Mixed second order inputs
        UiUj = zeros(size(i_UiUj,1),1);
        for i = 1:size(i_UiUj,1)
                UiUj(i,:) = u(i_UiUj(i,1),jj) .* u(i_UiUj(i,2),jj);
        end
    else
        i_UiUj = [];
    end
    n_uy = max(N,lengthu);
    
    mii=[];
    for i=1:min(N,lengthu)
        mii=[mii;[i i]];
    end
    iUY = nchoosek(1:n_uy,2);                            % Combination of indices if n(u)=n(measured outputs)=n_uy
    iUY=[iUY;flip(iUY,2)];
    [~,iUY1]=sort(iUY(:,1));
    iUY=[mii;iUY(iUY1,:)];
    if n_uy==N                                          % Eliminate indices aloted to u greater than n(u)
        elim = find(iUY(:,1)>lengthu);
        iUY(elim,:) = [];
    else                                                % Eliminate indices alotted to y greater than n(measured outputs)
        elim = find(iUY(:,2)>N);
        iUY(elim,:) = [];
    end
    UY = zeros(size(iUY,1),1);
    for i=1:size(iUY,1)
        UY(i,:) = u(iUY(i,1),jj).*XXQ(iUY(i,2),:);     % Compute U*Y
    end
    U2 = u(:,jj) .* u(:,jj);
    LogU = log(u(:,jj));                                % Natural log transformation
    ExpU = exp(-u(:,jj));                               % Exponential transformation
    InvU = 1./u(:,jj);                                  % Inverse of inputs
    InvUsq = (u(:,jj).*u(:,jj)).\1;                     % Inverse of squared inputs
    sqrtU = sqrt(u(:,jj));                              % Square root of inputs
    sqrtInvU = sqrtU.\1;                                % 1/sqrt(u)
    sigmoidU = (1+exp(-u(:,jj))).\1;                    % Sigmoid

    if lengthu>1
        Utrans = [ones(1,1); u(:,jj); UY; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
    else
        Utrans = [ones(1,1); u(:,jj); UY; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
    end
    Ntrans = size(Utrans,1);                            % NUMBER OF TRANSFORMATIONs
    
    % Arranging Utrans according to rank
    Utrans_r = Utrans(idOb,:);                          

    if Ntrans>85
        Ntrans = 85;
        Utrans_r = Utrans_r(1:Ntrans,:);
    end

    % Model parameters for particular rank
    sln = solnrank;
    indx = find(Obtrnd(sln,:));
    Cm = Cmat(sln,:);                 
    Am = Amat(sln,:);  
    Cm_ = reshape(Cm,Ntrans,N)';      
    A = reshape(Am,N,N)';
    C = Cm_(:,indx);      
    
    % Model simulation
    Utrans_ = Utrans_r(indx,:);                         % Selecting relevant input transformations
    if jj==1
        XQ(:,jj+1) = A*XXQ + C*Utrans_;
    else
        XQ(:,jj+1) = A*XQ(:,jj) + C*Utrans_;
    end
        
    XXQ = XQ(:,jj+1);
end

for j=1:size(XQ,1)
    y(j,:) = mnd(j) + 0.5*(XQ(j,:)-1)*(mxd(j)-mnd(j));
end

% Model interpretation
fprintf(['\nThe identified model is of the form \n \n y(k+1) = A*y(k) + C*Utrans(k) \n \n'])
fprintf(['Selected basis functions in Utrans for model ', num2str(sln), ': \n \n']);
BasisInterpretation(lengthu,iUY,i_UiUj,Obtrnd,idOb,sln,true);

fprintf(['\n Model parameters: \n']);
A
C


return

