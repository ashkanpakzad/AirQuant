function tblout = AQ_vertcat(t1, t2)
% for vertical concatenation of tables with different columns

%variable names
t1_names = t1.Properties.VariableNames;
t2_names = t2.Properties.VariableNames;

%shared variable names
t2_not_t1 = setdiff(t2_names, t1_names);
t1_not_t2 = setdiff(t1_names,t2_names);

% make cols in t1
for var = t2_not_t1
    t1.(var{1}) = cell(height(t1),1);
end

% make cols in t2
for var = t1_not_t2
    t2.(var{1}) = cell(height(t2),1);
end


%catenation
tblout = [t1; t2];

end