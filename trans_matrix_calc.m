function [P]=trans_matrix_calc(X,norder)

% %if n(i,1)==1
%         X1=importdata('activities/direct_care.txt');
% %elseif n(i,1)==2
%         X2=importdata('activities/housekeeping.txt');
% %elseif n(i,1)==3
%         X3=importdata('activities/mealtimes.txt');
% %elseif n(i,1)==4
%         X4=importdata('activities/medication_round.txt');
% %elseif n(i,1)==5
%         X5=importdata('activities/miscellaneous.txt');
% %elseif n(i,1)==6
%         X6=importdata('activities/personal.txt');
% %end


%X(isnan(X))=0;
%X=importdata('seq.dat');
N = max(max(X));                                   %# Number of states
[P, Q] = meshgrid(1:N, 1:N);
Y = [X, zeros(size(X, 1), 1)]';                    %# Pad for concatenation
if norder==1
    count_func = @(p, q)numel(strfind(Y(:)', [p, q])); %# Counts p->q transitions
else
    count_func = @(p, q)numel(strfind(Y(:)', [p, q, p]));
end
P = reshape(arrayfun(count_func, P, Q), N, N);
P = bsxfun(@rdivide, P, sum(P,2));

% N = max(max(X));                                     %# Number of states
% [p, q] = meshgrid(1:N, 1:N);
% Y = [X, zeros(size(X, 1), 1)]';                      %# Pad for concatenation
% count_func = @(u)numel(strfind(Y(:)', [p(u), q(u)])); %# Counts p->q transitions
% P = reshape(arrayfun(count_func, 1:numel(p)), N, N);

