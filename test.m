
% a = [1 3 5; 2 4 6; 7 8 10]

%for a = 10:20 
   %fprintf('value of a: %d\n', a);
%end

%for a = 1.0: -0.1: 0.0
   %disp(a)
%end

%v = [ 1; 2; 3; 4; 5; 6];	% creating a column vector of 6 elements
%v(3)


n = 1;
%area=zeros(1,n);

for i = 1:5
    x = rand;
    %keyboard;
    SN(i).b = i + 1;
    SN(i).a = x;
    SN(i).c = i + 10;
end

%SN(1).a
S = orderfields(SN);
S.a

%SNARRAY = struct2cell(SN);
%celldisp(SNARRAY)
%SNMAT = cell2mat(SNARRAY);

%OUT = sortrows(SNMAT,[1,3])



%T = struct2cell(SN); % convert the struct array to a table
%TM = cell2mat(T)

%sortedT = sortrows(T, 'b'); % sort the table by 'DOB'
%sortedS = table2struct(sortedT) % change it back to struct array if necessary




