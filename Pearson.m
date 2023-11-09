function [DataMatrix,DataName] = Pearson(Data)
DataMatrix = double(Data(2:end,:));
DataName = Data(1,:);
i = 1;
while i <= size(DataMatrix,2)
    PearsonList = zeros(1,size(DataMatrix,2));
    for j = i+1:size(DataMatrix,2)
         YData = DataMatrix(:,i);
         XData = DataMatrix(:,j);
         PearsonList(1,j)=corr(XData,YData,'type','Pearson');  
    end
    PLPosition = (abs(PearsonList) > 0.6);
    DataMatrix(:,PLPosition) = [];
    DataName(:,PLPosition) = [];
    i = i+1;
end