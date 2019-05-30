% Function to calculate the CFU values after care based on either
% uni-directional pick-up or bi-directional (type 2). To access, type
% either spread_cfu(sequence vector,model type)
%[spread,A,beta,lambda,V]=spread_cfu(seq,model,room,ACH,position,P) 
%P= position at which HCW is currently working
function [spread,A,beta,lambda,V]=spread_cfu(seq,model,room,ACH,position,P)   
%import=which particle room data to import 1 or 4
%ACH which column 4 or 6
if nargin==1
    model=1;
end
        [r,c]=size(seq);
        A=surface_area(r,c);
        %normrnd(18.6,4.1753,size(seq)); %Brouwer hand surface area
        lambda=lambda_transfer(r,c,seq);%normrnd(0.35,0.14^2,size(seq));%Transfer efficiency from surf to hand
        beta=beta_transfer(r,c,seq);%ph=zeros(m,1);
        V=zeros(size(seq));  %produces V values specific to data based 
 % cd G:\YAB\CFD For old cfucm2 data (wrong formula excel)
%   import=input('Which room?');
  if room==1
      v=importdata('cfucm2_yab.txt');%average cfu/cm2 values from CFD
      P=1;
  else
      v=importdata('cfucm2_hbn4.txt');
      %position=1;%input('Which Release point?');
  end
  
%    cd C:\Users\Phi\Documents\MATLAB
%   ACH=input('Which ACH rate?');
  if room==1 && ACH==4
      y=1;
  elseif room==1 && ACH==6
      y=2;
  elseif room==4 && ACH==4 && position==1
      y=1;
  elseif room==4 && ACH==4 && position==2
      y=3;
  elseif room==4 && ACH==4 && position==3
      y=5;
  elseif room==4 && ACH==4 && position==4
      y=7;
  elseif room==4 && ACH==6 && position==1
      y=2;
  elseif room==4 && ACH==6 && position==2
      y=4;
  elseif room==4 && ACH==6 && position==3
      y=6;
  elseif room==4 && ACH==6 && position==4
      y=8;
  end
  if room==4
      s=std(v(:,y)); %work out std from cfd because it doesn't produce one
  else
      s=std(v(:,2));
  end
        for i=1:length(seq)
            for j=1:size(seq,2)
                if seq(i,j)==1     %Equipment	std(v(:,y))
                    
                    MU = log(v(1*P,y)^2 / sqrt(s+v(1*P,y)^2));
                    SIGMA = sqrt(log(s/v(1*P,y)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);%normrnd(x(1,y),0.01,1,1);
                elseif seq(i,j)==2 %Patient	
                    MU = log(v(2*P,y)^2 / sqrt(s+v(2*P,y)^2));
                    SIGMA = sqrt(log(s/v(2,y)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);
                elseif seq(i,j)==3 %Hygiene products	
                    MU = log(v(3*P,y)^2 / sqrt(s+v(3*P,y)^2));
                    SIGMA = sqrt(log(s/v(3*P,y)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);
                elseif seq(i,j)==4 %Near bed objects	
                    MU = log(v(4*P,y)^2 / sqrt(s+v(4*P,y)^2));
                    SIGMA = sqrt(log(s/v(4*P,y)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);
                elseif seq(i,j)==5 %Far objects
                     MU = log(v(5*P,y)^2 / sqrt(s+v(5*P,y)^2));
                    SIGMA = sqrt(log(s/v(5*P,y)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);
                end
            end
        end
                
        %V=randi(5,size(seq)); %Dancer's threshold
       %V=round(rand(size(seq))); %<1
        %V=poissrnd(3,size(seq)); %surface pathogen loading cfu/cm^2
%   t = cputime;      
if model==1
    spread = ones(size(seq)) .* (  V.*lambda.*A.*(seq==1) +...
           V.*lambda.*A.*(seq==2)+ V.*lambda.*A.*(seq==3) + ...
           V.*lambda.*A.*(seq==4) + V.*lambda.*A.*(seq==5) + ...
           zeros(size(seq)).*( (seq==0) ));
else 

for i=1:length(seq)
     for j=1:size(seq,2)
        if seq(i,j)==1
            if j>=2  
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j));
            end
        elseif seq(i,j)==2
            if j>=2  
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j));
            end
        elseif seq(i,j)==3
            if j>=2  
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j));
            end

        elseif seq(i,j)==4
            if j>=2  
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j));
            end
        elseif seq(i,j)==5
            if j>=2  
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=(V(i,j)*lambda(i,j)*A(i,j));
            end

        elseif seq(i,j)==0
            spread(i,j)=0;
        end
    end
end

   
     %  spread(i,:)=spread(i,:)-rand*spread(;%dep(i,:,id);
   %current CFU=previous CFU - deposition+ previous accretion

end

       
spread=cumsum(spread,2);       
% e = cputime-t

 %{
 fig
 hold on
  plot(sum(seq~=0,2),spread(:,end),'bd');
  if room==1
      title(['YAB single room under','', num2str(ACH),'ACH']);
  else
      title(['HBN04 standard 4-bedded room under', num2str(ACH),'ACH']);
  end
  xlabel('Number of contacts n')
  ylabel('Colony forming units')
%}
%   b=boxplot([reshape(spread_YAB_4ACH,28000,1),reshape(spread_YAB_6ACH,28000,1)],'labels',{'4ACH','6ACH'});
%   ylabel('Colony forming units')
%   fig
%   boxplot([reshape(spread_yab_4,28000,1),reshape(spread_yab_6,28000,1)],'labels',{'4ACH','6ACH'})

%   boxplot([spread_1(:,end) spread_2(:,end)],'labels',{'No deposition','With deposition'});
%  fig
%  plot(sum(seq~=0,2),spread_1(:,end),'ok');
%  hold on
%   plot(sum(seq~=0,2),spread_2(:,end),'.r');
%   xlabel('Number of contacts n')%   ylabel('Colony forming units')

% fig
% subplot(1,2,1)
% [f x]=hist(reshape(spread_1,28000,1),1:40:400);
% bar(x,f/sum(f),'barwidth',0.5,'facecolor','r');
% xlabel('Colony forming units (Y)')
% ylabel('Frequency density')
% axis([-10 400 0 0.5])
% [f x]=hist(reshape(spread_2,28000,1),1:40:400);
% 
% subplot(1,2,2)
% bar(x,f/sum(f),'barwidth',0.5,'facecolor','r');
% xlabel('Colony forming units (Y)')
% ylabel('Frequency density')
% axis([-10 400 0 0.5])

