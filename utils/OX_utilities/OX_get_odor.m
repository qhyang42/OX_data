function [odorvec, category] = OX_get_odor(subjname) 
%%% get odor and category info from all available runs 
cuelistdir = '/Users/qhyang/Desktop/OX_DATA/cuelist';
odorvec = []; 
category = {}; 
for sesidx = 1: 16 % 16 sessions at most
    for runidx = 1:10
        try
            load(fullfile(cuelistdir, subjname, ['session', num2str(sesidx)], ['cuelist_sess', num2str(sesidx), '_run', num2str(runidx), '.mat']));
        catch
            continue
        end

    odorvec = [odorvec; cuelist.odor]; 
    category = [category, cuelist.category ]; 
    end 
end 

end 