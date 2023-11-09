%% program for spiking neural P systems-inspired evolutionary algorithm (SNPEA)
clear;clc;

% problem
%%pearson
DataBraintumor_ALL = importdata('braintumor.csv');
ClassifyNum = DataBraintumor_ALL(:,end);
DataBraintumor_ALL = DataBraintumor_ALL(:,1:(size(DataBraintumor_ALL,2)-1));
NameBraintumor = load('braintumor_gene_name.mat');
NameBraintumor = NameBraintumor.braintumor_gene_name;
Data = [NameBraintumor;DataBraintumor_ALL];
[DataMatrix,DataName] = Pearson(Data);
DataMatrix = [DataMatrix ClassifyNum];
save('DataMatrix.mat','DataMatrix','DataName');
gene_name = DataName;
TimeList = [];
AccuracyList = [];
MaxIteration = 300;
%%
for numofgene = 10:10
    for numofloop = 1:10
        tic
        maxfit = 0;
        % parameter setting
        Ps=50; % population size
        LengthOfProblem=length(gene_name);
        historialBestSolution=-9999*ones(1,Ps);
        historialBestBit=zeros(Ps,LengthOfProblem);
        globalmaxfit=-9999;
        
        
        % initialization
        P=rand(Ps,LengthOfProblem);

        
        %------------------------------------------------------------------%
        flagtermination=-99;
        Iteration=0;
        unchangedSolution=globalmaxfit;
        unchangedTime=0;
        globalmaxfitlist = [];
        globalmaxbitlist = [];
        globalmaxgenelist =[];
        eachgenmaxfit = [];
        eachgenlocation = [];
        eachgenname = [];
        fenjie = 25;
        svmlist = [];
        rflist = [];
        while Iteration < MaxIteration
            Iteration = Iteration+1;
            if Iteration == 1
                C1 = rand(Ps/2 ,LengthOfProblem);
                C = [C1 ; C1];
                [B] = pop_observe(C,numofgene);
                fitness=[];
                for i = 1:Ps
                    if i <= fenjie
                        acc = use_svm(B(i,1:LengthOfProblem),DataMatrix);
                        fitness = [fitness ; acc];
                    else
                        acc = use_knn(B(i , 1:LengthOfProblem),DataMatrix);
                        fitness = [fitness ; acc];
                    end
                end
            else
                [B] = pop_observe(P,numofgene);
                fitness=[];
                parfor i=1:Ps
                    if i <= fenjie
                        acc = use_svm(B(i,1:LengthOfProblem),DataMatrix);
                        fitness = [fitness ; acc];
                    else
                        acc = use_knn(B(i,1:LengthOfProblem),DataMatrix);
                        fitness = [fitness ; acc]; 
                    end
                end
            end
            
            %------------------------------------------------------%
            [maxfit,maxindex]=max(fitness);
            maxbit=B(maxindex,:);
            if maxindex <= fenjie
                if fenjie ==35
                    fenjie = fenjie;
                else
                    fenjie = fenjie + 1;  
                end
                svmlist = [svmlist ; Iteration];
            else
                if fenjie ==15
                    fenjie = fenjie;
                else
                    fenjie = fenjie - 1;   
                end
                rflist = [rflist ; Iteration];
            end
            eachgenmaxfit = [eachgenmaxfit ; maxfit];
            eachlocation = find(maxbit);
            eachgenlocation = [eachgenlocation ; eachlocation];
            eachgenge = gene_name(eachlocation);
            eachgenname = [eachgenname ; eachgenge];
            
            
            %-----------------------------------------------------%
            for i=1:Ps
                if historialBestSolution(1,i)<fitness(i,1)
                    historialBestSolution(1,i)=fitness(i,1);
                    historialBestBit(i,:)=B(i,:);
                end
            end
            %--------------------------------------------------%
            if (globalmaxfit<=maxfit)
                globalmaxfit=maxfit;
                globalmaxbit=maxbit;
                globalmaxfitlist = [globalmaxfitlist ; globalmaxfit];
                bestlocation = find(globalmaxbit);
                globalmaxbitlist = [globalmaxbitlist ; bestlocation];
                globalmaxgene = gene_name(bestlocation);
                globalmaxgenelist = [globalmaxgenelist ; globalmaxgene];
            else
                globalmaxfitlist = [globalmaxfitlist ; globalmaxfit];
                bestlocation = find(globalmaxbit);
                globalmaxbitlist = [globalmaxbitlist ; bestlocation];
                globalmaxgene = gene_name(bestlocation);
                globalmaxgenelist = [globalmaxgenelist ; globalmaxgene];
            end
            if maxfit == 100
                break
            end
            %-----------------------------------------------------%
            updateP=P;
            for i=1:Ps
                if unchangedTime<50
                    Plearn=0.05+0.15*rand;
                else
                    Plearn=0.2+0.8*rand;
                end
                for j=1:LengthOfProblem
                    if rand<Plearn
                        k1=ceil(Ps*rand);
                        while k1==i
                            k1=ceil(Ps*rand);
                        end
                        k2=ceil(Ps*rand);
                        while (k2==i)|(k2==k1)
                            k2=ceil(Ps*rand);
                        end
                        if historialBestSolution(1,k1)>historialBestSolution(1,k2)
                            b=historialBestBit(k1,j);
                        else
                            b=historialBestBit(k2,j);
                        end
                        if B(i,j)~=b

                            if (b>0.5)
                                updateP(i,j)=0.5 + 0.5*P(i,j);
                            else
                                updateP(i,j)=0.5*P(i,j);
                            end
                        end
                        clear b;
                    else
                        if B(i,j)~=globalmaxbit(1,j)
                            if globalmaxbit(1,j)>0.5
                                updateP(i,j)=0.5 + 0.5*P(i,j);
                            else
                                updateP(i,j)=0.5*P(i,j);
                            end
                        end
                    end
                end
            end
            P=updateP;
            
            bestprofit=globalmaxfit;
            if bestprofit>unchangedSolution
                unchangedSolution=bestprofit;
                unchangedTime=0;
            else
                unchangedTime=unchangedTime+1;
            end
            if unchangedTime>500000
                flagtermination=99;
            end

        end
 
        toc
        time = toc;
        save(['genedata_gene',num2str(numofgene-1),'_',num2str(numofloop),'.mat'],'globalmaxfitlist','globalmaxfit','time','Iteration')
        TimeList = [TimeList;time];
        AccuracyList = [AccuracyList;globalmaxfit];
        disp(['Time:',num2str(time),'*****Accuracy:',num2str(globalmaxfit)]);
    end
end
for i = 1:size(TimeList,1)
    disp(['Time:',num2str(TimeList(i,:)),'*****Accuracy:',num2str(AccuracyList(i,:))]);
    filename = 'data.txt';
    fid = fopen(filename,'a');
    fprintf(fid,'Time:%f*****Accuracy:%f\n',TimeList(i,:),AccuracyList(i,:));
end
[SizeR,~] = size(TimeList);
AverTime = sum(TimeList)/SizeR;
AverAcc = sum(AccuracyList)/SizeR;
disp(['AverTime:',num2str(AverTime),'*****AverAcc:',num2str(AverAcc)]);
fprintf(fid,'AverTime:%f*****AverAcc:%f',AverTime,AverAcc);
fclose(fid);