set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)


region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
cd 'graph'
if ~exist('tunnel_cor', 'dir')
    mkdir('tunnel_cor');
end
cd 'tunnel_cor'
Anova_name=[];
Anova_cor=[];
for fi=1:9
    corr=[];
    title=region;
    for regi = 1:4
        temp_corr=[];
        temp=tab.rspin{1, fi}(1+(regi-1)*5:regi*5);
        index=find(~cellfun(@isempty,temp));
        n=numel(index);
        for i=1:n
            for j=i+1:n
                [r,lags]=xcorr(temp{1,index(i)},temp{1,index(j)},'normalized');
                temp_corr=[temp_corr,max(r)];
            end
        end
        corr(regi)=mean(temp_corr);
        Anova_cor=[Anova_cor, temp_corr];
        title(regi)=strcat(title(regi)," n=",num2str(n));
        Anova_name=[Anova_name; repmat(categorical(region(regi)),numel(temp_corr),1)];
    end
    cat_name=categorical(title);

    cat_name=reordercats(cat_name,cellstr(title));
    figure
    bar(cat_name,corr);
    ylabel ("average correlations")
    ylim([0,1]);
    saveas(gcf,tab.Fol{1, fi},'png')
end
[~,~,stats1] = anova1(Anova_cor',Anova_name);
[ff_c,ff_means]=multcompare(stats1,0.05,'on','','s');
cd ..\..

%% FF
all_cor=[];
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];

cd 'graph'
if ~exist('tunnel_cor', 'dir')
    mkdir('tunnel_cor');
end
cd 'tunnel_cor'
Anova_cor=[];
Anova_name=[];
titl=region;
for regi = 1:4
    all_n(regi)=0;
    corr=[];
    for fi=1:9
        temp_corr=[];
        temp=tab.rspin{1, fi}(1+(regi-1)*5:regi*5);
        index=find(~cellfun(@isempty,temp));
        n=numel(index);
        c_n=0;
        for i=1:n
            if Index_fd{fi, (regi-1)*5+index(i)}=='ff'
                for j=i+1:n
                    if Index_fd{fi, (regi-1)*5+index(j)}=='ff'
                        [r,lags]=xcorr(temp{1,index(i)},temp{1,index(j)},'normalized');
                        temp_corr=[temp_corr,max(r)];

                    end

                end
                c_n=c_n+1;
            end
        end
        corr=[corr,temp_corr];
        Anova_cor=[Anova_cor, temp_corr];
        Anova_name=[Anova_name; repmat(categorical(region(regi)),numel(temp_corr),1)];
        all_n(regi)=all_n(regi)+c_n;
    end


    all_cor(regi)=mean(corr);
    std_err(regi)=std(corr)/sqrt(numel(corr));
    titl(regi)=strcat(titl(regi)," n=",num2str(all_n(regi)));
end
cat_name=categorical(titl);
cat_name=reordercats(cat_name,cellstr(titl));
bar(cat_name,all_cor);
hold on;
er = errorbar(cat_name,all_cor,std_err);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel ("average correlations")
ylim([0,1]);
title ("Feed Forward")
hold off
saveas(gcf,'All corr FF','png')
cd ../..
[~,~,stats1] = anova1(Anova_cor',Anova_name);
stats1.means(2)=[];
stats1.n(2)=[];
[ff_c,ff_means]=multcompare(stats1,0.05,'on','','s');
%% FB
all_cor=[];
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
Anova_cor=[];
Anova_name=[];
cd 'graph'
if ~exist('tunnel_cor', 'dir')
    mkdir('tunnel_cor');
end
cd 'tunnel_cor'
titl=region;
for regi = 1:4
    all_n(regi)=0;
    corr=[];
    for fi=1:9
        temp_corr=[];
        temp=tab.rspin{1, fi}(1+(regi-1)*5:regi*5);
        index=find(~cellfun(@isempty,temp));
        n=numel(index);
        c_n=0;
        for i=1:n
            if Index_fd{fi, (regi-1)*5+index(i)}=='fb'
                for j=i+1:n
                    if Index_fd{fi, (regi-1)*5+index(j)}=='fb'
                        [r,lags]=xcorr(temp{1,index(i)},temp{1,index(j)},'normalized');
                        temp_corr=[temp_corr,max(r)];

                    end
                    c_n=c_n+1;
                end
            end
        end
        corr=[corr,temp_corr];
        all_n(regi)=all_n(regi)+c_n;
        Anova_cor=[Anova_cor, temp_corr];
        Anova_name=[Anova_name; repmat(categorical(region(regi)),numel(temp_corr),1)];
    end

    all_cor(regi)=mean(corr);
    std_err(regi)=std(corr)/sqrt(numel(corr));
    titl(regi)=strcat(titl(regi)," n=",num2str(all_n(regi)));
end
cat_name=categorical(titl);
cat_name=reordercats(cat_name,cellstr(titl));
bar(cat_name,all_cor);
hold on;
er = errorbar(cat_name,all_cor,std_err);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel ("average correlations")
ylim([0,1]);
title("Feedback")
hold off
saveas(gcf,'All corr FB','png')
cd ../..
[~,~,stats1] = anova1(Anova_cor',Anova_name);
stats1.means(1)=[];
stats1.n(1)=[];
[fb_c,fb_means]=multcompare(stats1,0.05,'on','','s');
