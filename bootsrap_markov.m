% sequences = arrayfun(@(~) randi([1 5], [1 randi([500 1000])]), 1:400, ...
%     'UniformOutput',false)';
load('seq_yab.mat')
sequences=num2cell(seq_yab,2);
sequences = cellfun(@(x) nonzeros(x).', sequences, 'UniformOutput', false);
%# compute transition matrix from all sequences
trans = countFcn(sequences);

%# number of bootstrap samples to draw
Nboot = 1000;

%# estimate 95% confidence interval using bootstrapping
ci = bootci(Nboot, {@countFcn, sequences}, 'alpha',0.05);
ci = permute(ci, [2 3 1]);

% cil=sum(ci(:,:,1),2);
% ciu=sum(ci(:,:,2),2);
% 
% for i=1:5
% ci(i,:,1)=ci(i,:,1)./cil(i);   
% ci(i,:,2)=ci(i,:,2)./ciu(i);
% end

%# compute multiple transition matrices using bootstrapping
stat = bootstrp(Nboot, @countFcn, sequences);

%# display histogram for each entry in the transition matrix
sub = reshape(1:5*5,5,5);

for i=1:size(stat,2)
    subplot(5,5,sub(i))
    histogram(stat(:,i))
end

