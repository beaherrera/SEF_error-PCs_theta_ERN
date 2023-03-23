function data_mat = lock2event(data_cell, time_cell)
%% LOCK2EVENT epoch data [-.5 1] secs relative to event (time 0)

data_mat = cellfun(@(time, data) data(:, time >= -.5 & time <= 1.), ...
    time_cell, data_cell, "UniformOutput",false);
data_mat = cat(3, data_mat{:});

end