
clc;clear all;close all;
addpath('util');

load demo_accuracy5_BinaryOpt;

data_list_name = {'randn','randn','randn', ...
    'sector','sector','sector',...
    '20news','20news','20news',...
    'caltech','caltech','caltech',...
    'e2006','e2006','e2006',...
    'tdt2','tdt2','tdt2',...
    };

    
mm_list = [1024 512 256  ];
nn_list = [512 1024 2048 ];

mm_list = [mm_list mm_list mm_list mm_list mm_list mm_list mm_list];
nn_list = [nn_list nn_list nn_list nn_list nn_list nn_list nn_list];

for idat = 1:length(data_list_name)
    his1 = result{idat}.his1;
    his2 = result{idat}.his2;
    his3 = result{idat}.his3;
    his4 = result{idat}.his4;
    his5 = result{idat}.his5;
    all_avg = [Get2Decimal(his1(1));Get2Decimal(his2(1));Get2Decimal(his3(1));Get2Decimal(his4(1));Get2Decimal(his5(1))];
    all_std = [Get2Decimal(his1(2));Get2Decimal(his2(2));Get2Decimal(his3(2));Get2Decimal(his4(2));Get2Decimal(his5(2))];
    rank = RankList(all_avg);
    
    DataStr  = data_list_name{idat};
    fprintf('%s-%d-%d',DataStr,mm_list(idat),nn_list(idat));
    for iii=1:5
        if(rank(iii)==0)
            fprintf('& %.2f $\\pm$ %.2f',all_avg(iii),all_std(iii) );
        elseif(rank(iii)==1)
            fprintf('& \\cone{%.2f $\\pm$ %.2f }',all_avg(iii),all_std(iii));
        elseif(rank(iii)==2)
            fprintf('& \\ctwo{%.2f $\\pm$ %.2f }',all_avg(iii),all_std(iii));
        elseif(rank(iii)==3)
            fprintf('& \\cthree{%.2f $\\pm$ %.2f }',all_avg(iii),all_std(iii));
        end
    end
    fprintf('\\\\\n');
end

