function baselineinds = setbaseline(time, baseline_start, baseline_end)

temp_times = time-baseline_start;
[~,baselinetimestart_index] = min(abs(temp_times));
temp_times = time-baseline_end;
[~,baselinetimeend_index] = min(abs(temp_times));
%store in a single variable
baselineinds = baselinetimestart_index:baselinetimeend_index;
