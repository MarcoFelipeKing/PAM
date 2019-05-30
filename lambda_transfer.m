function lambda=lambda_transfer(r,c,seq)
lambda=zeros(r,c);
%randsample([40.0000,   41.0000,   27.6000,   38.4700,   33.4700,   65.8000],1,'true')
for i=1:r
    for j=1:c
        if seq(i,j)==2
            lambda(i,j)=0.35-rand*0.35;
        else
            lambda(i,j)=gamrnd(13.65,3)/100;
            
            if lambda(i,j)>1
               lambda(i,j)=1;
            end
        end
    end
end