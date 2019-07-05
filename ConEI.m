function ConMat=ConEI(Ne,Ni,pee,pii,pei,pie)
%ConMat=ConEI(800,200,0.05,0.3,0.3,0.3);

%Initialize interconnectivity matrix and intraconnectivity matrix
MatEE=zeros(Ne);
MatII=zeros(Ni);
MatEI=zeros(Ne,Ni);
MatIE=zeros(Ni,Ne);

%Generate (Ne+Ni) by (Ne+Ni) random matrix
%compare with corresponding probability to decide whether each pair of
%cells are connected
RandMat=rand(Ne+Ni);
for i=1:Ne
    for j=1:Ne
        if RandMat(i,j) <= pee
            MatEE(i,j)=1;
        end
    end
end

for i=1:Ne
    for j=1:Ni
        if RandMat(i,j+Ne) <= pei
            MatEI(i,j)=1;
        end
    end
end

for i=1:Ni
    for j=1:Ne
        if RandMat(i+Ne,j) <= pie
            MatIE(i,j)=1;
        end
    end
end

for i=1:Ni
    for j=1:Ni
        if RandMat(i+Ne,j+Ne) <= pii
            MatII(i,j)=1;
        end
    end
end

%Combine all connectivity matrices, 1 represents connected
ConMat=[MatEE,MatEI;MatIE,MatII];

%Set diagonal entries to 0 so that each cell doesn't connect with itself
for i=1:(Ne+Ni)
    ConMat(i,i)=0;
end

end
