% ManifoldData


for i=1:length(ManifoldData)-1

    SpkWindow_Nosepoke_All=vertcat(ManifoldData{i,2},ManifoldData{i+1,2});
    SpkWindow_Reward_All=vertcat(ManifoldData{i,3},ManifoldData{i+1,3});

end


