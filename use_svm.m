
function acc = use_svm(lie,DataMatrix)

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


[bestacc, bestc,bestg] = SVMcgForClass( t_train, p_train,-2,4,-4,4,5,0.5,0.5,0.9);
cmd = ['-t 0', '-c', num2str(bestc), '-g', num2str(bestg)];
model = libsvmtrain(t_train, p_train, cmd);

T_sim2 = libsvmpredict(t_test , p_test , model);


error2 = sum((T_sim2' == T_test )) / N * 100;
acc = error2 ;


