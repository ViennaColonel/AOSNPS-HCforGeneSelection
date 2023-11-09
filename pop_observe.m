function [B]=pop_observe(pop,numofgene)
A=pop;  
rowlengthofA=size(A,1);

B=(rand(size(pop))<A);
for i=1:rowlengthofA
    C = find(B(i,:)==1);
    rowofA=A(i,C);
    sortofrowofA=sort(rowofA,'descend');
    threshholdvalue=sortofrowofA(numofgene);
    referencevector=linspace(threshholdvalue,threshholdvalue,size(rowofA,2));
    D=(referencevector<rowofA);
    B(i,:) = 0;
    B(i,C) = D;
end