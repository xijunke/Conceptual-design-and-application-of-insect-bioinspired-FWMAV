%% optimal_aeroF_coeff
tic
% matlabpool open 4
% function obj_function=objfunction_aeroF_coeff(x)
ObjFunction=@objfunction_aeroF_stdev;  % fitness and constraint functions——调用目标函数objfunction_aeroF_coeff
nvars=3;   % Number of variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Call ga——遗传算法
% options=gaoptimset('Display','diagnose','Display','iter','StallTimeLimit',7200,'PopulationSize',100,...
%                'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
%                'UseParallel','always'); %并行计算—可行
% % options=gaoptimset('PopulationSize',100,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr},'UseParallel','always'); %并行计算可行 % 20150407
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fminsearchoptions = optimset('Display','iter','MaxFunEvals',1000,'MaxIter',1000,'TolCon',1e-6,'TolFun',1e-6,'UseParallel','always');
% % 'Display','notify-detailed'
% fminsearchoptions = optimset('Display','iter-detailed','MaxFunEvals',1000,'MaxIter',1000,'TolCon',1e-15,'TolFun',1e-15,'UseParallel','always');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options =gaoptimset(options,'Hybridfcn',{@fminsearch,fminsearchoptions});
% [x,fval,exitflag,output,population,scores]=ga(ObjFunction,nvars,[],[],[],[],[],[],[],options)
% matlabpool close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Call patternsearch——模式搜索
x0=[1,1,1];       % stdev =41.7446;
% psoptimset——Create pattern search options structure
% % Complete poll around current iterate; 'off'
options=psoptimset('CompletePoll','on','CompleteSearch','on','Display','iter','UseParallel','always'); 
[x,fval,exitflag,output] = patternsearch(ObjFunction,x0,[],[],[],[],[],[],[],options)
% matlabpool close
optimal_variablex=x;
object_value=fval;
save('result_x.txt', 'optimal_variablex', '-ASCII')
save('result_fval.txt', 'object_value', '-ASCII')
% output_1=[population,scores];
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\aerodynamic_moment\aerodynamic_moment_model2_2_FCoeff_optimal\population_scores.xlsx',output_1,'sheet1','A1:K100');
elapsedTime =toc/60;  
disp(['Caculation finished and the elapsedTime= ' num2str(elapsedTime ,'%5.4f') 'minutes'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x =[1.0000,1.0000,0.7500,0.5643,1.5000];  % fval =35.5047;  % Elapsed time is 583.787292 seconds.
% x =[1.0000,1.6003,0.7500,0.5574,0.0100];  % fval =35.5099;   % 26 minutes
% x= [1.0000,0.5000,0.3634, 0.5000];             %  fval =36.8167; % 7.3954minutes
% x =[1.0000,0.5000,0.3358,0.5000];              % fval=36.5981;  % 7.9058minutes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% output_1=[population,scores];
% xlswrite('D:\KXJ\hybrid_GA_fminsearch_WingM4_4\population_scores.xlsx',output_1,'sheet1','A1:J100');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


