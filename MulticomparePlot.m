CW09B.Properties.VariableNames{'Snowball'} = 'Porter';
CW09B.Properties.VariableNames{'SelSnowball'} = 'SelPorter';
groupNames = CW09B.Properties.VariableNames(:,2:end)
[p,tbl,stats] = friedman(table2array(CW09B(:,2:end)))
stats.gnames = groupNames;
[c,m,h,gnames] = multcompare(stats,'alpha',0.1,"ctype","hsd")

yticklabels(flip(groupNames))
title('CW09B')
ax = gca;
ax.FontSize = 14;
set(ax.Children,'LineWidth',2)
set(h,'Position',[0 0 800 800])