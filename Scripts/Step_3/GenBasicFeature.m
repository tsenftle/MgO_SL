function [Data, Head, Unit] = GenBasicFeature(Data_in, Head_in, Unit_in, Operation, Remain)

% GenBasicFeature is a script to generate secondary descriptors from
% primary descriptors by using the function "BasicOperation" and for-loop.

Data =[];
Head =[];
Unit =[];

for i = 1: length(Operation)
    [D_temp, H_temp, U_temp] = BasicOperation(Data_in, Head_in, Unit_in, char(Operation(i)));
    Data = [Data D_temp];
    Head = [Head H_temp];
    Unit = [Unit U_temp];
end

if Remain == 1
    Data = [Data_in Data];
    Head = [Head_in Head];
    Unit = [Unit_in Unit];
end

end