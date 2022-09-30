%% create ridgeplots
figure();set(gcf,'color',[1 1 1],'position',[0 0 width/3 height])
toi = find(cleaned.time{1}==.1):find(cleaned.time{1}==.25);
conds = {1;2;3;4;5;6;7};
load('vik.mat')
cols = {[0 0 1];...
        [0 0 .83];...
        [0 0 .66];...
        [.33 0 .33];...
        [.66 0 0];...
        [.83 0 0];...
        [1 0 0]};
% cols = {colourFriendlyMaps('blue');...
%         colourFriendlyMaps('blue');...
%         colourFriendlyMaps('blue');...
%         colourFriendlyMaps('black');...
%         colourFriendlyMaps('reddishpurple');...
%         colourFriendlyMaps('reddishpurple');...
%         colourFriendlyMaps('reddishpurple')};
labels = { {'\itp\rm=0'},...
           {'\itp\rm=.17'},...
           {'\itp\rm=.33'},...
           {'\itp\rm=.5'},...
           {'\itp\rm=.67'},...
           {'\itp\rm=.83'},...
           {'\itp\rm=1'}};
xlims=[-1 4];
for i = 1:7
    
    data = squeeze(mean(mean(mean(alldata(:,conds{i},el,toi),3),4),2));
    
    [f,xi]=ksdensity(data);
    h=subplot(7,1,i);
    hold on;
    p=fill(xi,f,[.9 .9 .9]);
    p.FaceColor = cols{i};
    p.EdgeColor = [.3 .3 .3];
    p.FaceAlpha = .3;
    p.LineWidth = 2;
    xlim(xlims)
    box off
    plot(xlims,[0 0],'color',[.3 .3 .3],'linewidth',2)
    ax=gca;
    ax.YColor='none';
    ax.Color = [1 1 1];
    if i ~= 7
        ax.XColor = 'none';
    else
        ax.XColor = [.3 .3 .3];
    end
    ax.FontSize=16;
    t=text(xlims(1)+.1,.5,1,labels{i},'fontsize',16)
    t.Interpreter='tex';
%     set(h,'Position',[.1 .9 .4 .4])
%     axis off;
end
%%
