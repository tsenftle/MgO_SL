function [Data, Head, Unit] = GenSelfFeature(Data_in, Head_in, Unit_in, Operation, Remain)

% GenSelfFeature is a script to generate secondary descriptors from
% primary descriptors by using the function "SelfOperation" and for-loop.

Data =[];
Head =[];
Unit =[];

for i = 1: length(Operation)
    for j = 1: size(Data_in, 2)
        [D_temp, H_temp, U_temp] = SelfOperation(Data_in(:,j), Head_in(j), Unit_in(j), char(Operation(i)));
        Data = [Data D_temp];
        Head = [Head H_temp];
        Unit = [Unit U_temp];
    end
end

if Remain == 1
    Data = [Data_in Data];
    Head = [Head_in Head];
    Unit = [Unit_in Unit];
end

end