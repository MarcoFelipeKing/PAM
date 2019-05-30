%# function test_spread_cfu just for calculating spread without room based on v.dat 
%# [spread,A,lambda,beta,V]=test_spread_cfu(seq,model)
%# New model: YiAi/Ah=lambdai*Vi*Ai/As+(1-beta)*Yi_1*Ai/Ah
%#
function [spread,A,lambda,beta,V]=new_spread_cfu(seq,model)

if nargin==1
    model=2;
end
spread=zeros(size(seq));
[r,c]=size(seq);
        A=surface_area(r,c);
        %normrnd(18.6,4.1753,size(seq)); %Brouwer hand surface area
        lambda=lambda_transfer(r,c,seq);%normrnd(0.35,0.14^2,size(seq));%Transfer efficiency from surf to hand
        beta=beta_transfer(r,c,seq);%ph=zeros(m,1);
        %v=importdata('v.dat')';    % Test room data
        %v=importdata('activities/v_1_exp.txt')';%cfu from single room
        %v=importdata('activities/v_N2_2_exp.txt')';
        
        %in old model V was cfu/cm^2
        %Lets assume they are total cfu on surface to get units correct
        vv=importdata('cfucm2_yab.txt');% YAB data 
        v=zeros(size(vv)); % for 4 on v(:,1) and 6 ach v(:,2)
        v(:,1)=vv(:,2);
        v(:,2)=std(vv(:,1));
        
        %v=importdata('V_single.txt');
        %v=importdata('V_N2.txt');
        %v=importdata('V_N1.txt');
        V=zeros(size(seq));  %produces V values specific to data based 
        %cd G:\YAB\CFD
        %v=importdata('cfucm2_yab.txt');
        for i=1:length(seq)
            for j=1:size(seq,2)
                if seq(i,j)==1     %Equipment	
                    MU = log(v(1,1)^2 / sqrt(v(1,2)+v(1,1)^2));
                    SIGMA = sqrt(log(v(1,2)/v(1,1)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);%normrnd(x(1,y),0.01,1,1);%poissrnd(4,1);
                elseif seq(i,j)==2 %Patient	
                    MU = log(v(2,1)^2 / sqrt(v(2,2)+v(2,1)^2));
                    SIGMA = sqrt(log(v(2,2)/v(2,1)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);%normrnd(x(2,y),0.01,1,1);%poissrnd(6,1);
                elseif seq(i,j)==3 %Hygiene products	
                    MU = log(v(3,1)^2 / sqrt(v(3,2)+v(3,1)^2));
                    SIGMA = sqrt(log(v(3,2)/v(3,1)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);%normrnd(x(3,y),0.01,1,1);%poissrnd(7,1);
                elseif seq(i,j)==4 %Near bed objects	
                    MU = log(v(4,1)^2 / sqrt(v(4,2)+v(4,1)^2));
                    SIGMA = sqrt(log(v(4,2)/v(4,1)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);%=normrnd(x(4,y),0.01,1,1);%poissrnd(9,1);
                elseif seq(i,j)==5 %Far objects
                    MU = log(v(5,1)^2 / sqrt(v(5,2)+v(5,1)^2));
                    SIGMA = sqrt(log(v(5,2)/v(5,1)^2 + 1));
                    V(i,j)=lognrnd(MU,SIGMA);%normrnd(x(5,y),0.01,1,1);%poissrnd(6,1);
                else
                    V(i,j)=0;
                end
            end
        end
        
if model==1 %This doesn't assume any differences in V A or lambda. Ignores beta
    spread = ones(size(seq)) .* (  V.*lambda.*A.*(seq==1) +...
           V.*lambda.*A.*(seq==2)+ V.*lambda.*A.*(seq==3) + ...
           V.*lambda.*A.*(seq==4) + V.*lambda.*A.*(seq==5) + ...
           zeros(size(seq)).*( (seq==0) ));
else 
            A_H=54; %Total hand size. https://journals.sagepub.com/doi/abs/10.1177/154193128603000417
for i=1:length(seq)
     for j=1:size(seq,2)
        if seq(i,j)==1
            A_S=10; %Total surface size
            if j>=2  
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S+(1-beta(i,j))*spread(i,j-1);
            %(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S;
            end
        elseif seq(i,j)==2
            A_S=10; %Total surface size
            if j>=2  
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S+(1-beta(i,j))*spread(i,j-1);
            %(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=%(V(i,j)*lambda(i,j)*A(i,j));
            end
        elseif seq(i,j)==3
            A_S=10; %Total surface size
            if j>=2  
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S+(1-beta(i,j))*spread(i,j-1);
            %(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S;%(V(i,j)*lambda(i,j)*A(i,j));
            end

        elseif seq(i,j)==4
            A_S=10; %Total surface size
            if j>=2  
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S+(1-beta(i,j))*spread(i,j-1);
            %(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S;%(V(i,j)*lambda(i,j)*A(i,j));
            end
        elseif seq(i,j)==5
            A_S=10; %Total surface size
            if j>=2  
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S+(1-beta(i,j))*spread(i,j-1);
            %(V(i,j)*lambda(i,j)*A(i,j))-beta(i,j)*spread(i,j-1);  
            else
            spread(i,j)=V(i,j)*lambda(i,j)*A_H/A_S;%(V(i,j)*lambda(i,j)*A(i,j));
            end

        elseif seq(i,j)==0
            spread(i,j)=0;
        end
    end
end
spread=cumsum(spread,2);  
end