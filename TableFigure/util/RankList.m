function [rank] = RankList(fs)
n = length(fs);
[fs_sort] = sort(fs);
fs_sort = unique(fs_sort);
if(length(fs_sort)==1)
    rank = ones(n,1);
elseif(length(fs_sort)==2)
    min_f_1 = fs_sort(1);
    min_f_2 = fs_sort(2);
    min_f_3 = fs_sort(2);
    rank = zeros(n,1);
    rank(find(fs==min_f_1))=1;
    rank(find(fs==min_f_2))=2;
    rank(find(fs==min_f_3))=2;
elseif(length(fs_sort)>=3)
    min_f_1 = fs_sort(1);
    min_f_2 = fs_sort(2);
    min_f_3 = fs_sort(3);
    rank = zeros(n,1);
    rank(find(fs==min_f_1))=1;
    rank(find(fs==min_f_2))=2;
    rank(find(fs==min_f_3))=3;
    
end
