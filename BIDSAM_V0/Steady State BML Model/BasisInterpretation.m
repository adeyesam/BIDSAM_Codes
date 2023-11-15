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
        % fprintf('Selected basis functions \n');
        for i=1:numel(idsol)
            if idsol(i)==1
                fprintf('1\n');
            elseif idsol(i)>Ut_i(1) && idsol(i)<=Ut_i(2)
                fprintf(['u', num2str(idsol(i)-Ut_i(1)), '\n'])    
            elseif idsol(i)>Ut_i(2) && idsol(i)<=Ut_i(3)
                fprintf(['u', num2str(i_UiUj(idsol(i)-Ut_i(2),1)), '*u', num2str(i_UiUj(idsol(i)-Ut_i(2),2)), '\n'])    
            elseif idsol(i)>Ut_i(3) && idsol(i)<=Ut_i(4)
                fprintf(['u', num2str(idsol(i)-Ut_i(3)), '^2 \n'])   
            elseif idsol(i)>Ut_i(4) && idsol(i)<=Ut_i(5)
                fprintf(['log(u', num2str(idsol(i)-Ut_i(4)), ')\n'])  
            elseif idsol(i)>Ut_i(5) && idsol(i)<=Ut_i(6)
                fprintf(['exp(-u', num2str(idsol(i)-Ut_i(5)), ') \n'])  
            elseif idsol(i)>Ut_i(6) && idsol(i)<=Ut_i(7)
                fprintf(['1/u', num2str(idsol(i)-Ut_i(6)), '\n'])    
            elseif idsol(i)>Ut_i(7) && idsol(i)<=Ut_i(8)
                fprintf(['1/u', num2str(idsol(i)-Ut_i(7)), '^2 \n'])   
            elseif idsol(i)>Ut_i(8) && idsol(i)<=Ut_i(9)
                fprintf(['sqrt(u', num2str(idsol(i)-Ut_i(8)), ') \n'])    
            elseif idsol(i)>Ut_i(9) && idsol(i)<=Ut_i(10)
                fprintf(['1/sqrt(u', num2str(idsol(i)-Ut_i(9)), ') \n']) 
            elseif idsol(i)>Ut_i(10) && idsol(i)<=Ut_i(11)
                fprintf(['1/(1+exp(-u', num2str(idsol(i)-Ut_i(10)), ')) \n'])   
            end
        end
    end
end