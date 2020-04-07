function [Data, Head, Unit] = SelfOperation(Data_in, Head_in, Unit_in, Operation)

% SelfOperation is a script to generate secondary descriptors from
% primary descriptors by using mathematical operations, such as: 
% ^-1, ^0.5, ^2 and ^3, and other operation, such as: exp, log and abs.

% Operations: {'^r', '^I', '^2', '^3', 'log', 'exp', 'abs'}
Data = zeros(size(Data_in,1),1);
Head = cell(1,1);
Unit = zeros(1,1);

switch Operation
    case '^r'
        if sum(Data_in)==0
            Data = Data_in;
        else
            Data = abs(Data_in).^0.5;
        end
        Head = strcat('(',Head_in,')^0.5');
        Unit = Unit_in/2;
    case '^I'
        if sum(Data_in)==0
            Data = Data_in;
        else
            Data = Data_in.^-1;
        end
        Head = strcat('(',Head_in,')^-1');
        Unit = -Unit_in;
    case '^2'
        if sum(Data_in)==0
            Data = Data_in;
        else
            Data = Data_in.^2;
        end
        Head = strcat('(',Head_in,')^2');
        Unit = Unit_in*2;
    case '^3'
        if sum(Data_in)==0
            Data = Data_in;
        else
            Data = Data_in.^3;
        end
        Head = strcat('(',Head_in,')^3');
        Unit = Unit_in*3;
    case 'log'
        if Unit_in == 0
            if sum(Data_in)==0
                Data = Data_in;
            else
                Data = log(abs(Data_in));
            end
            Head = strcat('log(',Head_in,')');
            Unit = NaN;
        end
    case 'exp'
        if Unit_in == 0
            if sum(Data_in)==0
                Data = Data_in;
            else
                Data = log(abs(Data_in));
            end
            Head = strcat('log(',Head_in,')');
            Unit = NaN;
        end
    case 'abs'
        if sum(Data_in)==0
            Data = Data_in;
        else
            Data = abs(Data_in);
        end
        Head = strcat('|',Head_in,'|');
        Unit = Unit_in;
    otherwise
        disp('unknown operation')
end
% Remove zero matrix
unique= find(any(Data));
Data = Data(:, unique);
Head = Head(:, unique);
Unit = Unit(:, unique);

end