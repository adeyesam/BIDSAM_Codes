function [] = BasisInterpretation(lengthu,i_UiUj,Obtrnd,idOb,sln,flag11)

    Ut_i(1) = 1;                                % Constant is 1:Ut_i(1)
    Ut_i(2) = Ut_i(1) + lengthu;                % u is Ut_i(1)+1:Ut_i(2)
    Ut_i(3) = Ut_i(2) + size(i_UiUj,1);         % ui*uj is Ut_i(3)+1:Ut_i(4)
    Ut_i(4) = Ut_i(3) +lengthu;                 % u^2 is Ut_i(4)+1:Ut_i(5)
    Ut_i(5) = Ut_i(4) +lengthu;                 % log(u) is Ut_i(5)+1:Ut_i(6)
    Ut_i(6) = Ut_i(5) +lengthu;                 % exp(-u) is Ut_i(6)+1:Ut_i(4)
    Ut_i(7) = Ut_i(6) +lengthu;                 % inv(u) is Ut_i(7)+1:Ut_i(8)
    Ut_i(8) = Ut_i(7) +lengthu;                 % inv(u^2) is Ut_i(8)+1:Ut_i(9)
    Ut_i(9) = Ut_i(8) +lengthu;                 % sqrt(u) is Ut_i(9)+1:Ut_i(10)
    Ut_i(10) = Ut_i(9) +lengthu;                % 1/sqrt(u) is Ut_i10)+1:Ut_i(11)
    Ut_i(11) = Ut_i(10) +lengthu;               % sigmoid(1/(1-exp(-u)) is Ut_i(11)+1:Ut_i(12)
    
    if flag11
        id_sol = find(Obtrnd(sln,:));
        idsol = idOb(id_sol);
        fprintf('Basis function           Indices \n');
        for i=1:numel(idsol)
            if idsol(i)==1
                Const=1;
                fprintf('Constant                     %d \n',1);
            elseif idsol(i)>Ut_i(1) && idsol(i)<=Ut_i(2)
                ui=idsol(i)-Ut_i(1);
                fprintf('ui                          %d \n',idsol(i)-Ut_i(1));
            elseif idsol(i)>Ut_i(2) && idsol(i)<=Ut_i(3)
                uiuj=i_UiUj(idsol(i)-Ut_i(2),:);
                fprintf('uiuj                      %d   %d \n',i_UiUj(idsol(i)-Ut_i(2),:));
            elseif idsol(i)>Ut_i(3) && idsol(i)<=Ut_i(4)
                u_sq=idsol(i)-Ut_i(3);
                fprintf('u_sq                        %d \n',idsol(i)-Ut_i(3));
            elseif idsol(i)>Ut_i(4) && idsol(i)<=Ut_i(5)
                logu=idsol(i)-Ut_i(4);
                fprintf('log(u)                      %d \n',idsol(i)-Ut_i(4));
            elseif idsol(i)>Ut_i(5) && idsol(i)<=Ut_i(6)
                expu=idsol(i)-Ut_i(5);
                fprintf('exp(-u)                     %d \n',idsol(i)-Ut_i(5));
            elseif idsol(i)>Ut_i(6) && idsol(i)<=Ut_i(7)
                Invu=idsol(i)-Ut_i(6);
                fprintf('Inv(u)                     %d \n',idsol(i)-Ut_i(6));
            elseif idsol(i)>Ut_i(7) && idsol(i)<=Ut_i(8)
                Invu2=idsol(i)-Ut_i(7);
                fprintf('Inv(u^2)                     %d \n',idsol(i)-Ut_i(7));
            elseif idsol(i)>Ut_i(8) && idsol(i)<=Ut_i(9)
                sqrtu=idsol(i)-Ut_i(8);
                fprintf('sqrt(u)                     %d \n',idsol(i)-Ut_i(8));
            elseif idsol(i)>Ut_i(9) && idsol(i)<=Ut_i(10)
                Invsqrtu=idsol(i)-Ut_i(9);
                fprintf('1/sqrt(u)                     %d \n',idsol(i)-Ut_i(9));
            elseif idsol(i)>Ut_i(10) && idsol(i)<=Ut_i(11)
                sigmoidu=idsol(i)-Ut_i(10);
                fprintf('sigmoid(-u)                     %d \n',idsol(i)-Ut_i(10));
            end
        end
    end
end