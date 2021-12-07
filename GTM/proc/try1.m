t_len = length(prob_tight);
a = struct;
for ind = 1:t_len
    a(ind).q = sum(prob_tight(ind).res);
end
plot([a.q])


%% Functions: