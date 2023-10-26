function [Utranss, iUY, i_UiUj] = InputTransformation(data,uu,IUY,IUU,flag101)
% For 'first time', all indices of UY and UU to be calculated, afterwardss, no need to repeat calculations
    if flag101
        N = size(data,1);                                  
        lengthu = size(uu,1);                             
        lengthu2 = size(uu,2);                            
        lenreadpone = size(data,2);
    
        if lengthu>1
            i_UiUj = nchoosek(1:lengthu,2);                 % Indices for Mixed second order inputs
            UiUj = zeros(size(i_UiUj,1),lengthu2);
            for i = 1:size(i_UiUj,1)
                    UiUj(i,:) = uu(i_UiUj(i,1),:) .* uu(i_UiUj(i,2),:);
            end
        else
            i_UiUj = [];
        end
        
        n_uy = max(N,lengthu);
        mii=[];
        for i=1:min(N,lengthu)
            mii=[mii;[i i]];
        end
        iUY = nchoosek(1:n_uy,2);                           % Combination of indices if n(u)=n(measured outputs)=n_uy
        iUY=[iUY;flip(iUY,2)];
        [~,iUY1]=sort(iUY(:,1));
        iUY=[mii;iUY(iUY1,:)];
        if n_uy==N                                          % Eliminate indices aloted to u greater than n(u)
            elim = iUY(:,1)>lengthu;
            iUY(elim,:) = [];
        else                                                % Eliminate indices alotted to y greater than n(measured outputs)
            elim = iUY(:,2)>N;
            iUY(elim,:) = [];
        end
        UY = zeros(size(iUY,1),lenreadpone);
        for i=1:size(iUY,1)
            UY(i,:) = uu(iUY(i,1),:).*data(iUY(i,2),:);      % Compute U*Y
        end
    
        U2 = uu .* uu;
        LogU = log(uu);                                      % Natural log transformation
        ExpU = exp(-uu);                                     % Exponential transformation
        InvU = 1./uu;                                        % Inverse of inputs
        InvUsq = (uu.*uu).\1;                                % Inverse of squared inputs
        sqrtU = sqrt(uu);                                    % Square root of inputs
        sqrtInvU = sqrtU.\1;                                 % 1/sqrt(u)
        sigmoidU = (1+exp(-uu)).\1;                          % Sigmoid

        if lengthu>1
            Utranss = [ones(1,lengthu2); uu; UY; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        else
            Utranss = [ones(1,lengthu2); uu; UY; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        end
    else

        lengthu_2 = size(uu,2);
        iUY = IUY;
        i_UiUj = IUU;
        if ~isempty(i_UiUj)
            for i = 1:size(i_UiUj,1)
                    UiUj(i,:) = uu(i_UiUj(i,1),:) .* uu(i_UiUj(i,2),:);
            end
        end
        
        for i=1:size(iUY,1)
            UY(i,:) = uu(iUY(i,1),:).*data(iUY(i,2),:);      % Compute U*Y
        end
    
        U2 = uu .* uu;
        LogU = log(uu);                                      % Natural log transformation
        ExpU = exp(-uu);                                     % Exponential transformation
        InvU = 1./uu;                                        % Inverse of inputs
        InvUsq = (uu.*uu).\1;                                % Inverse of squared inputs
        sqrtU = sqrt(uu);                                    % Square root of inputs
        sqrtInvU = sqrtU.\1;                                 % 1/sqrt(u)
        sigmoidU = (1+exp(-uu)).\1;                          % Sigmoid
        
        if ~isempty(i_UiUj)
            Utranss = [ones(1,lengthu_2); uu; UY; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        else
            Utranss = [ones(1,lengthu_2); uu; UY; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        end

    end
end
