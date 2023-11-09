function acc = use_knn(lie,DataMatrix)
res = DataMatrix;
index = find(lie);

P_train = res([1: 7,11:17,21:27,31:33,35:38], index)';
T_train = res([1: 7,11:17,21:27,31:33,35:38], end)';
M = size(P_train, 2);

P_test = res(1: end, index)';
T_test = res(1: end, end)';
N = size(P_test, 2); 


[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input );
t_train = T_train;
t_test  = T_test ;


p_train = p_train'; p_test = p_test';
t_train = t_train'; t_test = t_test';
acclist =[];
for k = 3:2:15
    mdl = ClassificationKNN.fit(p_train,t_train,'NumNeighbors',k);
    label = predict(mdl,p_test);

    error2 = sum((label == t_test )) / N * 100;
    disp(['取近邻数K = ' num2str(k),'; 此时的准确率为 ' num2str(error2) '%'])
    acclist = [acclist , error2];
end
acc = max(acclist);
disp(acc)