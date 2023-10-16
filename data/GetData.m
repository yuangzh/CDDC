function [A,DataStr] = GetData(iwhich,type)

rand('seed',0); randn('seed',0);

if(type==1)% In type A, there is no restriction on m and n
    ms = [2048;1024;1024];
    ns = [1024;1024;2048];
else % In type B, there is a restriction on m and n that m>n
    ms = [2048;2048;2048];
    ns = [256;512;1024];
end


switch iwhich
    case 11 % nonsparse
        m = ms(1); n = ns(1);A = GetRand(m,n); DataStr = sprintf('randn-%d-%d',m,n);
    case 12
        m = ms(2); n = ns(2); A = GetRand(m,n); DataStr = sprintf('randn-%d-%d',m,n);
    case 13
        m = ms(3); n = ns(3); A = GetRand(m,n); DataStr = sprintf('randn-%d-%d',m,n);
        
    case 21 % nonsparse
        m = ms(1); n = ns(1);A = GetGisette(m,n); DataStr =  sprintf('Gisette-%d-%d',m,n);
    case 22
        m = ms(2); n = ns(2); A = GetGisette(m,n); DataStr =  sprintf('Gisette-%d-%d',m,n);
    case 23
        m = ms(3); n = ns(3); A = GetGisette(m,n); DataStr =  sprintf('Gisette-%d-%d',m,n);
        
    case 31 % sparse
        m = ms(1); n = ns(1);A = GetSector(m,n); DataStr = sprintf('sector-%d-%d',m,n);
    case 32
        m = ms(2); n = ns(2); A = GetSector(m,n); DataStr = sprintf('sector-%d-%d',m,n);
    case 33
        m = ms(3); n = ns(3); A = GetSector(m,n); DataStr = sprintf('sector-%d-%d',m,n);
        
    case 41 % sparse
        m = ms(1); n = ns(1); A = Get20Newsgroups(m,n); DataStr = sprintf('20news-%d-%d',m,n);
    case 42
        m = ms(2); n = ns(2); A = Get20Newsgroups(m,n); DataStr = sprintf('20news-%d-%d',m,n);
    case 43
        m = ms(3); n = ns(3); A = Get20Newsgroups(m,n); DataStr = sprintf('20news-%d-%d',m,n);
        
    case 51 % sparse
        m = ms(1); n = ns(1); A = GetRealDataE2006(m,n);  DataStr = sprintf('e2006-%d-%d',m,n);
    case 52
        m = ms(2); n = ns(2); A = GetRealDataE2006(m,n);  DataStr = sprintf('e2006-%d-%d',m,n);
    case 53
        m = ms(3); n = ns(3); A = GetRealDataE2006(m,n);  DataStr = sprintf('e2006-%d-%d',m,n);
        
    case 61 % sparse
        m = ms(1); n = ns(1); A = GetTDT2(m,n); DataStr =  sprintf('tdt2-%d-%d',m,n);
    case 62
        m = ms(2); n = ns(2); A = GetTDT2(m,n); DataStr =  sprintf('tdt2-%d-%d',m,n);
    case 63
        m = ms(3); n = ns(3); A = GetTDT2(m,n); DataStr =  sprintf('tdt2-%d-%d',m,n);
        
    case 71 % nonsparse
        m = ms(1); n = ns(1); A = GetCnn_4096d_Caltech(m,n); DataStr =  sprintf('CnnCaltech-%d-%d',m,n);
    case 72
        m = ms(2); n = ns(2); A = GetCnn_4096d_Caltech(m,n); DataStr =  sprintf('CnnCaltech-%d-%d',m,n);
    case 73
        m = ms(3); n = ns(3); A = GetCnn_4096d_Caltech(m,n); DataStr =  sprintf('CnnCaltech-%d-%d',m,n);
        
end
A = NormData(A,1); A = full(A);




function A = GetTDT2(m,n)
load TDT2
A = fea; A = double(A);
A = RandomSelRowsDelZeroCol(A,m);
A = RandomSelCol(A,n);
A = NormData(A,2);


function [A] = RandomSelRowsDelZeroCol(A,m)
[data_m,data_n] = size(A);
s1 = randperm(data_m,m);
A = A(s1,:);
zeroColumns = all(A == 0, 1);
A(:, zeroColumns) = [];

function [A] = RandomSelCol(A,n)
[data_m,data_n] = size(A);
s2 = randperm(data_n,n);
A = A(:,s2);

function A = GetSector(m,n)
load sector_train;
A = x;A = double(A);
A = RandomSelRowsDelZeroCol(A,m);
A = RandomSelCol(A,n);





function A = GetGisette(m,n)
load gisette;
A = x; A = double(A);
A = RandomSelRowsDelZeroCol(A,m);
A = RandomSelCol(A,n);
A = NormData(A,2);


function A = GetCifar(m,n)
load cifar;
A = data; A = double(A);
A = A / norm(A(:),inf);
A = RandomSelRowsDelZeroCol(A,m);
A = RandomSelCol(A,n);
A = NormData(A,1);

function A = GetCnn_4096d_Caltech(m,n)
load cnn_4096d_Caltech;
A = x; A = double(A);
A = RandomSelRowsDelZeroCol(A,m);
A = RandomSelCol(A,n);
A = NormData(A,2);


function A = GetRand(m,n)
A = randn(m,n);


function A = Get20Newsgroups(m,n)
load 20Newsgroups;
A = fea;A = double(A);
A = RandomSelRowsDelZeroCol(A,m);
A = RandomSelCol(A,n);




function A = GetRealDataE2006(m,n)
load E2006_5000_10000;
A = x;A = double(A); %[A] = center(A);
A = RandomSelRowsDelZeroCol(A,m);
A = RandomSelCol(A,n);
A = NormData(A,2);



function [A] = center(A)
[m,n] = size(A);
e = ones(m,1);
A = A - e*e'*A;



function [A] = removeZero(A)
% Find the row and column indices containing all zeros
zero_rows = all(A==0,2);
zero_cols = all(A==0,1);
A = A(~zero_rows, ~zero_cols);




function [A] = scaleA(A)
[m,n] = size(A);
seq = randperm(m*n,round(0.2*m*n));
A(seq) = A(seq)*100;
