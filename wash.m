function [cfu,ph]=wash(spread,seq,ACTid,r)
%# r is for parametric testing of increased handwashing
% [cfu,ph]=wash(spread,seq,ACTid,r)
cfu=ones(length(spread),1);
seqpat=sum(seq'==2); %by patient contact

%Prob Direct Care	Housekeeping	Mealtimes	Medication Round	Miscellaneous	Personal Care
wash_yab=[0.401015228	0.235294118	0.47619048	0.513513514	0.361111111	0.692307692];

%Prob type Handwash only	Gloves only	 Alcohol rub only
% prob_type=[0.258883249,	0.015228426,	0.152284264;
% 0.176470588,	0.235294118,	0.058823529;
% 0.19047619,	0,	0.333333333;
% 0.207207207,	0.027027027,	0.351351351;
% 0.222222222,	0,	0.180555556;
% 0.615384615,	0.153846154,	0.153846154;
%0.24361949,	0.027842227,	0.213457077;
%];%probability of hand hygiene based on patient and care type

prob_type=[0.607142857,	0.035714286,	0.357142857;
0.375,	0.5,	0.125;
0.363636364,	0,	0.636363636;
0.353846154,	0.046153846,	0.6;
0.551724138,	0,	0.448275862;
0.666666667,	0.166666667,	0.166666667;
];
seqlength=sum(seq'~=0); %by contact length


ph=zeros(length(spread),1);
randr=rand(length(spread),1);


    for i=1:length(spread)
    
        if numel(ACTid)==1
            ACTid=repmat(ACTid,length(spread),1);
        end
        d=ACTid;

        w=wash_yab(1,d);   
             
      
        if randr(i)<= w+r
            w_type=randsample(1:3,1,'true',prob_type(d(i),:)/sum(sum(prob_type(d(i),:),2)));
            if w_type==1,
                h=(1-10.^(-abs(normrnd(1.6175,0.1206)))); %handwash with water soap from log reductions sickbert bennet (don't 
            elseif w_type==3
                h=(1-10.^(-abs(normrnd(1.10,0.8129))));   %alc
            else
                h=1-gamrnd(5.91,0.40)*0.01;%Gloves Montville 2001  %gloves
            end
        else 
            h=0;
        end
        %h=randsample([h rand],1,'true',[0.8 0.2]);
        cfu(i,1)=spread(i,end)*(1-h);
        ph(i,1)=h;
        
    
    end
   