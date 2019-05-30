function lambda=lambda_transfer(r,c,seq)
lambda=zeros(r,c);
%randsample([40.0000,   41.0000,   27.6000,   38.4700,   33.4700,   65.8000],1,'true')
for i=1:r
    for j=1:c
        if seq(i,j)==2 %This must be because of surface category
            lambda(i,j)=normrnd(0.03,0.01); %Porous surfaces
        else
            lambda(i,j)=normrnd(0.69,0.13);%Non-porous surface From our paper %gamrnd(13.65,3)/100; 
            
            if lambda(i,j)>1
               lambda(i,j)=1;
            end
        end
    end
end