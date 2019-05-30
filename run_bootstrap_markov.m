 
for i=1:6
    if i==1
             X1=importdata('activities/direct_care.txt');
P1hat=trans_matrix_calc(X1,1);
Ptilde=P1hat./(sum(reshape(P1hat,25,1))+P1hat);
    elseif i==2
            X2=importdata('activities/housekeeping.txt');
P2hat=trans_matrix_calc(X2,1);
Ptilde=P2hat./(sum(reshape(P2hat,25,1))+P2hat);
    elseif i==3
        X3=importdata('activities/mealtimes.txt');
P3hat=trans_matrix_calc(X3,1);
Ptilde=P3hat./(sum(reshape(P3hat,25,1))+P3hat);
    elseif i==4
               X4=importdata('activities/medication_round.txt');
P4hat=trans_matrix_calc(X4,1);
Ptilde=P4hat./(sum(reshape(P4hat,25,1))+P4hat);
    elseif i==5
         X5=importdata('activities/miscellaneous.txt');
P5hat=trans_matrix_calc(X5,1);
Ptilde=P5hat./(sum(reshape(P5hat,25,1))+P5hat);
    elseif i==6
           X6=importdata('activities/personal.txt');
P6hat=trans_matrix_calc(X6,1);
Ptilde=P6hat./(sum(reshape(P6hat,25,1))+P6hat);

    end
PTILDE(:,:,i)=Ptilde; %write to array

[ci stat trans]=bootstrap_markov(Ptilde,500); %bootstrap on Ptilde, use trans
ci_cell{i}=ci;
stat_cell{i}=stat;
trans_cell{i}=trans;
Ptrans(:,:,i)=trans;

emis = ones(5,1);
for j=1:500
[~, seq_t(j,:)]=hmmgenerate(30,trans,emis);
[~, seq_L(j,:)]=hmmgenerate(30,ci(:,:,1),emis); %generates chain based on ci_L
[~, seq_U(j,:)]=hmmgenerate(30,ci(:,:,2),emis);
end

[spread(:,:,i),A,lambda,beta,V]=test_spread_cfu(seq_t,2);
% [spread_L,A,lambda,beta,V]=test_spread_cfu(seq_L,2);
% [spread_U,A,lambda,beta,V]=test_spread_cfu(seq_U,2);
% 
% STATS(i,1)=mean(abs(spread(:,end)-spread_L(:,end))./spread(:,end)*100);
% STATS(i,2)=std(abs(spread(:,end)-spread_L(:,end))./spread(:,end)*100);
% STATS(i,3)=mean(abs(spread(:,end)-spread_U(:,end))./spread(:,end)*100);
% STATS(i,4)=std(abs(spread(:,end)-spread_U(:,end))./spread(:,end)*100);

end




      

 



