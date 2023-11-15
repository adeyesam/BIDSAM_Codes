%% The file loads and simulates hierarchically ranked models (based on AICc+2Cve values)
% from the branch and bound algorithm. 'solnrank' is the rank of model we desire to simulate.
% The rows of 'inputs' are the input variables while the columns are the time instants


function y = simulate_model_steady(inputs, solnrank)

load('solution_workspace')
u = inputs;
[lengthu, lengthu2] = size(u);

% Scaling inputs 
for i = 1:lengthu
    u(i,:) = 1+2*(u(i,:)-umin(i))/(umax(i)-umin(i));
end


for jj=1:lengthu2
    
    % Compute input transformations
    i_UiUj = nchoosek(1:lengthu,2);                     % Indices for Mixed second order inputs
    UiUj = zeros(size(i_UiUj,1),1);
    for i = 1:size(i_UiUj,1)
            UiUj(i,:) = u(i_UiUj(i,1),jj) .* u(i_UiUj(i,2),jj);
    end
    U2 = u(:,jj) .* u(:,jj);
    LogU = log(u(:,jj));                                % Natural log transformation
    ExpU = exp(-u(:,jj));                               % Exponential transformation
    InvU = 1./u(:,jj);                                  % Inverse of inputs
    InvUsq = (u(:,jj).*u(:,jj)).\1;                     % Inverse of squared inputs
    sqrtU = sqrt(u(:,jj));                              % Square root of inputs
    sqrtInvU = sqrtU.\1;                                % 1/sqrt(u)
    sigmoidU = (1+exp(-u(:,jj))).\1;                    % Sigmoid

    Utrans = [1; u(:,jj); UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
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
    Cm=Cmat(sln,:);                 
    Cm_=reshape(Cm,Ntrans,N)';    
    C=Cm_(:,indx);      
    
    % Model simulation
    Utrans_ = Utrans_r(indx,:);                         % Selecting relevant input transformations
    yy = C*Utrans_;
    for j=1:size(yy,1)
        y(j,jj) = mnd(j) + 0.5*(yy(j,:)-1)*(mxd(j)-mnd(j));
    end
end

% Model interpretation
fprintf('\nThe identified model is of the form \n \n y = C*Utrans \n \n')
fprintf(['Selected basis functions in Utrans for model ', num2str(sln), ': \n \n']);
BasisInterpretation(lengthu,i_UiUj,Obtrnd,idOb,sln,true);

fprintf(['\n Model parameters: \n']);
C

return

