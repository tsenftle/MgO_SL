function [Data, Head, Unit] = RmVarZero(Data_in, Head_in, Unit_in, Tolerance)

for i = 1: size(Data_in,2)
    if var(Data_in(:,i)) <= Tolerance
        Data_in(:,i) = 0;
    end
end

unique= find(any(Data_in));
Data = Data_in(:, unique);
Head = Head_in(:, unique);
Unit = Unit_in(:, unique);

end