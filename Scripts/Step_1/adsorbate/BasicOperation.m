function [Data, Head, Unit] = BasicOperation(Data_in, Head_in, Unit_in, Operation)

% BasicOperation is a script to generate secondary descriptors from
% primary descriptors by using basic mathematical operations (plus, minus, 
% times and devision) and absolute operation.

% Operations: {'+', '-', '*', '/', '+abs', '-abs', '/abs', '*abs'}

index = 1: size(Data_in, 2);    % Number of descriptors in Data_in matrix
n_c=nchoosek(index,2);          % Number of possible combinations

Data = zeros(size(Data_in,1), size(n_c, 1));
Head = cell(1, size(n_c, 1));
Unit = zeros(1, size(n_c, 1));

for i = 1: size(n_c, 1)
    switch Operation
        case '+'
            if Unit_in(n_c(i,1)) == Unit_in(n_c(i,2))
                if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
                else
                    Data(:,i) = Data_in(:,n_c(i,1)) + Data_in(:,n_c(i,2));
                end
                Head(i) = strcat('(',Head_in(n_c(i,1)),'+',Head_in(n_c(i,2)),')');
                Unit(i) = Unit_in(n_c(i,1));
            end
        case '-'
            if Unit_in(n_c(i,1)) == Unit_in(n_c(i,2))
                if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
                else
                    Data(:,i) = Data_in(:,n_c(i,1)) - Data_in(:,n_c(i,2));
                end
                Head(i) = strcat('(',Head_in(n_c(i,1)),'-',Head_in(n_c(i,2)),')');
                Unit(i) = Unit_in(n_c(i,1));
            end
        case '*'
            if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
            else
                Data(:,i) = Data_in(:,n_c(i,1)) .* Data_in(:,n_c(i,2));
            end
            Head(i) = strcat('(',Head_in(n_c(i,1)),'*',Head_in(n_c(i,2)),')');
            Unit(i) = Unit_in(n_c(i,1)) + Unit_in(n_c(i,2));
        case '/'
            if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
            else
                Data(:,i) = Data_in(:,n_c(i,1)) ./ Data_in(:,n_c(i,2));
            end
            Head(i) = strcat('(',Head_in(n_c(i,1)),'/',Head_in(n_c(i,2)),')');
            Unit(i) = Unit_in(n_c(i,1)) - Unit_in(n_c(i,2));
        case '+abs'
            if Unit_in(n_c(i,1)) == Unit_in(n_c(i,2))
                if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
                else
                    Data(:,i) = abs(Data_in(:,n_c(i,1)) + Data_in(:,n_c(i,2)));
                end
                Head(i) = strcat('|',Head_in(n_c(i,1)),'+',Head_in(n_c(i,2)),'|');
                Unit(i) = Unit_in(n_c(i,1));
            end
        case '-abs'
            if Unit_in(n_c(i,1)) == Unit_in(n_c(i,2))
                if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
                else
                    Data(:,i) = abs(Data_in(:,n_c(i,1)) - Data_in(:,n_c(i,2)));
                end
                Head(i) = strcat('|',Head_in(n_c(i,1)),'-',Head_in(n_c(i,2)),'|');
                Unit(i) = Unit_in(n_c(i,1));
            end
        case '*abs'
            if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
            else
                Data(:,i) = abs(Data_in(:,n_c(i,1)) .* Data_in(:,n_c(i,2)));
            end
            Head(i) = strcat('|',Head_in(n_c(i,1)),'*',Head_in(n_c(i,2)),'|');
            Unit(i) = Unit_in(n_c(i,1)) + Unit_in(n_c(i,2));
        case '/abs'
            if (sum(Data_in(:,n_c(i,1)))==0) | (sum(Data_in(:,n_c(i,2)))==0)
                    Data(:,i) = Data(:,i);
            else
                Data(:,i) = abs(Data_in(:,n_c(i,1)) ./ Data_in(:,n_c(i,2)));
            end
            Head(i) = strcat('|',Head_in(n_c(i,1)),'/',Head_in(n_c(i,2)),'|');
            Unit(i) = Unit_in(n_c(i,1)) + Unit_in(n_c(i,2));
        otherwise
            disp('unknown operation')
    end
end

% Remove zero matrix
unique= find(any(Data));
Data = Data(:, unique);
Head = Head(:, unique);
Unit = Unit(:, unique);

end