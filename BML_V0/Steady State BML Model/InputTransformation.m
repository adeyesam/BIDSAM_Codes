function [Utranss, i_UiUj] = InputTransformation(data,uu,IUU,flag101)
% For 'first time', all indices of UY and UU to be calculated, afterwardss, no need to repeat calculations
    if flag101
        lengthu = size(uu,1);                           
        lengthu2 = size(uu,2);                             
    
        if lengthu>1
            i_UiUj = nchoosek(1:lengthu,2);                     % Indices for Mixed second order inputs
            UiUj = zeros(size(i_UiUj,1),lengthu2);
            for i = 1:size(i_UiUj,1)
                    UiUj(i,:) = uu(i_UiUj(i,1),:) .* uu(i_UiUj(i,2),:);
            end
        else
            i_UiUj = [];
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
            Utranss = [ones(1,lengthu2); uu; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
            % Utranss = [uu; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        else
            Utranss = [ones(1,lengthu2); uu; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
            % Utranss = [uu; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        end

    else

        lengthu_2 = size(uu,2);
        i_UiUj = IUU;
        for i = 1:size(i_UiUj,1)
                UiUj(i,:) = uu(i_UiUj(i,1),:) .* uu(i_UiUj(i,2),:);
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
            Utranss = [ones(1,lengthu_2); uu; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
           %  Utranss = [uu; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        else
            Utranss = [ones(1,lengthu_2); uu; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
            % Utranss = [uu; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
        end
        
    end
end
