clc;
clear;
%____________________________________________
%Input Parameters
%____________________________________________
D=7; %Interleaver Depth. More the value of D more is the randomness
n=1; %Number of Data blocks that you want to send.
data=[1 1 1 1 0 0 0 1 1 1 0 1 0 1 1 0 1 0 1 0 1 0 1 1]; % A constant pattern used as a data
%____________________________________________
%Make the length of the data to be a multiple of D. This is for
%demonstration only.Otherwise the empty spaces has to be zero filled.
data=repmat(data,1,n); %send n blocks of specified data pattern
if mod(length(data),D) ~=0
    data=[data data(1:(D*(fix(length(data)/D+1))-length(data)))];
end

%We are sending D blocks of similar data
%intlvrInput=repmat(data(1:n),[1 D]);
%fprintf('Input Data to the Interleaver -> \n');
%disp(char(intlvrInput));
%disp('____________________________________________________________________________');

%INTERLEAVER
%Writing into the interleaver row-by-row
permuterIndex=randperm(D);
intrlvrOutput=[];
index=1;

for i=1:fix(length(data)/D)
    intrlvrOutput=[intrlvrOutput data(permuterIndex+(i-1)*D)];
end

for i=1:mod(length(data),D)
    intrlvrOutput=[intrlvrOutput data(permuterIndex(i)+fix(length(data)/D)*D)];
end

uncorruptedIntrlvrOutput=intrlvrOutput;

%Corrupting the interleaver output by inserting 10 '*'
%intrlvrOutput(1,25:34)=zeros(1,10)+42;

%DEINTERLEAVER
deintrlvrOutput=[];

for i=1:fix(length(intrlvrOutput)/D)
    deintrlvrOutput(permuterIndex+(i-1)*D)=intrlvrOutput((i-1)*D+1:i*D);
end

for i=1:mod(length(intrlvrOutput),D)
    deintrlvrOutput((fix(length(intrlvrOutput)/D))*permuterIndex+i)=intrlvrOutput((i+1):(i+1)*D);
end

deintrlvrOutput=deintrlvrOutput;
disp('Given Data -->');
disp(data);
disp(' ')
disp('Permuter Index-->')
disp(permuterIndex);
disp(' ')
disp('PseudoRandom Interleaver Output -->');
disp(uncorruptedIntrlvrOutput);
disp(' ')
disp('Interleaver Output after being corrupted by 10 symbols of burst error - marked by ''*''->');
disp(char(intrlvrOutput));
disp(' ')
disp('PseudoRandom Deinterleaver Output -->');
disp(deintrlvrOutput);