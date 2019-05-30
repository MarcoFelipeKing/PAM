%%%%%%%%%%
% This is a PAM of student office at the University of Hong Kong
% Considers only right hand
% 29 students

%To do 23/05
% 1. Create histograms of number of contacts for each student
% 2. Generate Markov chains
% 3. Compare above with what is observed. How do you do this?
% 4. Calculate accretion using new PAM
%

%%%%
m=1e6; %Samples
studNum=1;%29; %Number of students Use 1 for General Matrix
max_contact_num=691;  %Maximum number of contacts between students arriving and leaving
emis = ones(34,1); %# get a sample of length = 5
seqXX=zeros(m,max_contact_num,studNum); %This defines the maximum number of contacts and students
moves=zeros(m,studNum); %Eeach student will do X moves
P=importdata('P_General.mat'); %imports matrices for each student P(:,:,k) k students
hist_pdf=importdata('hist_pdf_office_General.mat'); %histograms of each student's moves
hist_pdf(:,2)=[]; %removes
%% Pre-amble to make P out of frequencies in P_Transition_Fequencies.mat file

% for i=1:29
% P(:,:,i)=P_Transition_Fequencies(:,:,i)./sum(P_Transition_Fequencies(:,:,i),2);
% end
% P_Transition_Fequencies(isnan(P_Transition_Fequencies))=0;
% P=P_Transition_Fequencies;

%%%
%THis makes the general hist_pdf, but you need a specific one for each
%student
%h=histogram(eachtimecount(:,2),0:1:700);
%hist_office(:,1)=h.BinEdges;
%hist_office(:,2)=h.Values;






for i=1:1
%define max number of contacts a student makes each time
moves(:,i)=randsample(1:length(hist_pdf),m,'true',hist_pdf(:,2));%randsample(max_contact_num,m,'true',hist_pdf(:,:,i)./sum(hist_pdf(:,:,i),2));
end

%% Runs student tracking
for n=1:studNum % For each student
% X1=importdata('activities/direct_care.txt');
% P1=trans_matrix_calc(X1,1); %this makes the transision matrix from touch
% data in direct_care.txt
% P=P1;


% picks the value of moves for the Markov chain to run because there is no absorbing state.
%moves(:,n)=randsample(1:25,m,'true',hist_pdf{n}./sum(hist_pdf{n}))';
% hist(moves,0:2:25)



for i=1:m %replicas of each student
    
    [~,track] = hmmgenerate(moves(i,1), P(:,:,n), emis);%hmmgenerate(moves(n,1), P(:,:,n), emis);
    %if track = 1
    %seq(i,1:moves(i,n))=track;
      seqXX(i,1:moves(i,n),n)=track;
end


%seqXX(:,:,n)=seq;

end

seqXX(seqXX==34)=0; %change position of leave

for i=1:10
    stairs(seqXX(i,:))
hold on
end
axis([0 200 0 34])


%Plot the mode
plot(mode(seqXX(1:10000,:)))
axis([0 50 0 34])

%stairs(mean(seqXX),'k-','LineWidth',2)
%plot(hist_pdf(:,1),hist_pdf(:,2),'o')
%histogram(seqXX) 
%Find non zero entries per row
figure
histogram(sum(seqXX~=0,2),'normalization','probability','BinEdges',1:700);
hold on
histogram(eachtimecount(:,2),'normalization','probability','BinEdges',1:700);%'NumBins',50);

axis([0 1000 0 0.2])
hold off

mean(sum(seqXX~=0,2))
mean(eachtimecount(:,2))

(inv(eye(33)-P(1:33,1:33)))*ones(33,1)

%% Student office
for n=1:studNum %number of students
[Y,~,~,~,~]=new_spread_cfu(seqXX(:,:,n),2);
end

%% This is OLD YAB
%Single room summation 1..4
 
 seqXX=importdata('seqXX.mat');
 for ii=1:6 %Number of care types
     [Y,~,~,~,~]=new_spread_cfu(seqXX(:,:,ii),2);
     for pat=1:4 % Number of patients
         [Y1(:,pat,ii),~]=wash(Y(:,end),seqXX(:,:,ii),ii,0);
          Y(:,end)=Y1(:,pat,ii);
     end
     
 end

 %%
 %Multi-bed room summation 1..4
 for pos=1:4 %position of infectious patient

     for ii=1:6 %category of care carried out
            %[Y,~,~,~,~]=spread_cfu(seqXX(:,:,ii),2,4,6,pos,pat);%first pat
         
    
    for pat=1:4
         if pat==1
             [Y,~,~,~,~]=spread_cfu(seqXX(:,:,ii),2,4,6,pos,pat); %Pat referes to the patient at which the HCW is working
         [Y4(:,pat,ii,pos),~]=wash(Y(:,end),seqXX(:,:,ii),ii,0);
         else
             [Y,~,~,~,~]=spread_cfu(seqXX(:,:,ii),2,4,6,pos,pat); %Pat referes to the patient at which the HCW is working
         [Y4(:,pat,ii,pos),~]=wash(Y(:,end)+Y4(:,pat-1,ii,pos),seqXX(:,:,ii),ii,0);
             
         end
    end
    
     end
     
 end
    
% %[seqX1,~]=transition_matrix_seq(1,1);
%  for i=1:6
%      for j=4:2:6  % comparing 4 vs 6 in single room
%  [Y,~,~,~,~]=test_spread_cfu(seqXX(:,:,1),2);
%  [Y1_6,h]=wash(Y(:,end),seqX1,i,0);
%      end
%  end
%  YW(:,1)=Y1;
%  YW(:,2)=Y1_6;
%   
%  YW(:,3)=Y2;
%  YW(:,4)=Y2_6;
%  
%  YW(:,5)=Y3;
%  YW(:,6)=Y3_6;
%  
%  YW(:,7)=Y4;
%  YW(:,8)=Y4_6;
%  
%  YW(:,9)=Y5;
%  YW(:,10)=Y5_6;
%  
%  YW(:,11)=Y6;
%  YW(:,12)=Y6_6;
 
  %%
 % comparing 4 ACH vs 6 in single room
 for j=4:2:6
     for i=1:6
         [Y,~,~,~,~]=spread_cfu(seqXX(:,:,i),2,1,j,1,1);%first pat
         [Yw,~]=wash(Y(:,end),seqXX(:,:,i),i,0);
         if j==4
             YW1_46(:,2*i-1)=Yw(:,end);
         else
             YW1_46(:,2*i)=Yw(:,end);
         end
     end
     
 end
 
  boxplot(YW1_46/mean(YW1_46(:,1)), {reshape(repmat('A':'F',2,1),12,1) repmat((1:2)',6,1)} ,'factorgap',10,'color','rk','symbol','r')
    set(gca,'ytick',0:5)
    set(gca,'yticklabel',0:5)
    set(gca,'xtick',1.5:3.1:19)
    axis([0 19 -1 5])
    
set(gca,'xticklabel',{'Direct care','Housekeeping','Mealtimes','Medication','Miscellaneous','Personal care'})
ylabel('Normalised Y');
legend(findobj(gca,'Tag','Box'),{'6 ach^{-1}','4 ach^{-1}'})
 %%
 %for i=1:4
     [Y,~,~,~,~]=spread_cfu(seqX1,2,4,6,1,1);%first pat
     [Yw,h]=wash(Y(:,end),seqX1,1,0);
     YW4(:,1)=Yw;
     
     [Y,~,~,~,~]=spread_cfu(seqX1,2,4,6,i,1); %second pat
     [Yw+Y,h]=wash(Y(:,end),seqX1,1,0);
 
 
 figure('color',[1 1 1]) % Single room Y1
 set(0,'DefaultAxesFontName', 'Times New Roman')
for i=1:6
    subplot(3,2,i)
    boxplot(Y1(:,:,i)./mean(Y1(:,1,1)),'color','r')
    axis([0  5 -1 10])
    set(gca,'ytick',0:2:10)
    %set(gca,'yticklabel',0:2:10)
    %set(gca,'xticklabel',{'1','2','3','4'})
ylabel('Normalised Y');
end


%%
 figure('color',[1 1 1]) %plot multibed
 set(0,'DefaultAxesFontName', 'Times New Roman')
for i=1:6
    subplot(3,2,i)
    boxplot(Y4(:,:,i,1)./mean(Y1(:,1,1)),'color','k')
    axis([0  5 -1 10])
    set(gca,'ytick',0:2:10)
    %set(gca,'yticklabel',0:2:10)
    %set(gca,'xticklabel',{'1','2','3','4'})
ylabel('Normalised Y');
end
h=gcf;
set(h,'PaperOrientation','landscape');

print(gcf, '-dpdf', 'Y4_1.pdf');

figure('color',[1 1 1]) %plot multibed averaged
 set(0,'DefaultAxesFontName', 'Times New Roman')
for i=1:6
    subplot(3,2,i)
    boxplot(Y4_mean(:,:,i)./mean(Y1(:,1,1)),'color','k')
    axis([0  5 -1 10])
    set(gca,'ytick',0:2:10)
    %set(gca,'yticklabel',0:2:10)
    %set(gca,'xticklabel',{'1','2','3','4'})
ylabel('Normalised Y');
end
h=gcf;
set(h,'PaperOrientation','portrait');

print(gcf, '-dpdf', 'Y4_Y11.pdf');


%%
for i=1:4 %moving Y1 and Y4 to boxplottable matrix
    YWW(:,2*i-1,:)=Y1(:,i,:);
    YWW(:,2*i,:)=Y4(:,i,:,1);
end
    
figure('color',[1 1 1]) %plot multibed averaged side by side
 set(0,'DefaultAxesFontName', 'Times New Roman')
 for i=1:6
     subplot(3,2,i)
     boxplot(YWW(:,:,i)./mean(YWW(:,1,1)), {reshape(repmat('A':'D',2,1),8,1) repmat((1:2)',4,1)} ,'factorgap',10,'color','rk','symbol','')
ylabel('Normalised Y');
%legend(findobj(gca,'Tag','Box'),'Four-bed','Single-bed')
 set(gca,'ytick',0:9)
    set(gca,'yticklabel',0:2:10)
    set(gca,'xtick',1.5:2.7:10)
    axis([0 12 -1 5])
    set(gca,'xticklabel',{'Patient 1','Patient 2','Patient 3','Patient 4'});
 end    


