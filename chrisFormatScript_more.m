m=1;
lEvt = cell(length(data.LL_TS),1);
for i = 1:length(data.LL_TS)
    if data.LL_Rf_Mask(i) == 0
        lEvt{m} = 'LL_U';
        m = m+1;
    end
    if data.LL_Rf_Mask(i) == 1
        lEvt{m} = 'LL_R';
        m = m+1;
    end
end

m=1;
rEvt = cell(length(data.RL_TS),1);
for i = 1:length(data.RL_TS)
    if data.RL_Rf_Mask(i) == 0
        rEvt{m} = 'RL_U';
        m = m+1;
    end
    if data.RL_Rf_Mask(i) == 1
        rEvt{m} = 'LL_R';
        m = m+1;
    end
end

Evt = [lEvt; rEvt];
EvtTm = [data.LL_TS; data.RL_TS];
[l,k] = sort(EvtTm);
Evt = Evt(k);