function beta=beta_transfer(r,c,seq)
beta=zeros(r,c);
hs=importdata('hs.dat'); %Montville data for hand-to-spigot
%randsample([40.0000,   41.0000,   27.6000,   38.4700,   33.4700,   65.8000],1,'true')
for i=1:r
    for j=1:c
        if seq(i,j)==2
            beta(i,j)=normrnd(0.8,0.13);% Porous0.35-rand*0.35; %this is rusin's data
        else
            beta(i,j)=10.^randsample(hs(:,2),1,'true')/100;%Non-porous - Montville's data
            
            if beta(i,j)>1
               beta(i,j)=1;
            end
        end
    end
end