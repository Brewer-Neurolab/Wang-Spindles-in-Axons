%% Length FF
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
regLabel1ff=[];regLabel1fb=[];
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
well_regList={"EC-DG","DG-CA3","CA3-CA1","CA1-EC"};
figOrder=[1 2 4 3];
if ~exist('graph', 'dir')
    mkdir('graph');
end

cd 'graph'
R_square=[];
r=[];

moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','mode','sd','n'});
figure( 'Position', [100 100 700 600])
binEdge = logspace(2,3.5,50); binCenter = convert_edges_2_centers(binEdge);
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];
    cross_anova_ff=[];
    c=0;
    n=0;
    for fi=1:9
        temp_setff=[];temp_setfb=[];
        tmp=tab.Len{1, fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff'
                temp_setff=[temp_setff,tmp{k+(regi-1)*5}.*0.04];
                if ~isempty(tmp{k+(regi-1)*5})
                    c=c+1;
                end
            elseif Index_fd{fi,k+(regi-1)*5}=='fb'
                temp_setfb=[temp_setfb,tmp{k+(regi-1)*5}.*0.04];

            end
        end
        totaln=numel(temp_setfb)+numel(temp_setff);
        [hist_objectff,hist_plotff,histcountff(coliff,:)] = log_binned_histogram(temp_setff, binEdge, "pdf",0,totaln);
        coliff=coliff+1;
        n=n+numel(temp_setff);
        r=[r,temp_setff];
        cross_anova_ff=[cross_anova_ff,temp_setff];
        regLabel1ff = [regLabel1ff; repmat(categorical(region(regi)),numel(temp_setff),1)];
    end

    histcountff(~any(~isnan(histcountff), 2),:)=[];
    meanCff = mean(histcountff,1); stdCff = stdErr(histcountff,1);
    meanCff=meanCff(~isnan(meanCff));
    binCenterff=binCenter(~isnan(meanCff));

    if regi~=2
        Length_ff(regi)={cross_anova_ff};
    else
        Length_ff(regi)={[]};
    end


    if regi~=2 & regi~=4
        [f_ff,gof_ff]=fit(log(binCenterff'), meanCff', 'gauss1');
        mu_ff(regi) = f_ff.b1; sigma_ff(regi) = f_ff.c1;
        [~,a2]=max(f_ff(log(binCenterff)));
        disp('mod')
        disp(binCenterff(a2))
        t_val=f_ff(log(binCenterff))'.*binCenterff;
        f_mean=mean(t_val);
        [~,a2]=min(abs(t_val-f_mean));
        disp('mean')
        disp(binCenterff(a2))
        disp('median')
        disp(exp(mu_ff(regi)))
        disp('std')
        disp(exp(sigma_ff(regi)))

        disp(c);disp(n);
    elseif regi~=2
        [f_ff,gof_ff]=fit(log(binCenterff'), meanCff', 'gauss2');
        mu_ff1 = f_ff.b1;
        mu_ff2 = f_ff.b2;
        sigma_ff1= f_ff.c1/sqrt(2);
        sigma_ff2= f_ff.c2/sqrt(2);
        disp('median')
        disp(exp(mu_ff1))
        disp(exp(mu_ff2))
        disp('std')
        disp(exp(sigma_ff1))
        disp(exp(sigma_ff2))
        disp(c);disp(n);
    else 
        [f_ff,gof_ff]=fit(log(binCenterff'), meanCff', 'gauss2','lower',[0,0,0,0,6.4,0],'upper',[Inf,Inf,Inf,Inf,7.2,Inf]);
        mu_ff1 = f_ff.b1;
        mu_ff2 = f_ff.b2;
        sigma_ff1= f_ff.c1/sqrt(2);
        sigma_ff2= f_ff.c2/sqrt(2);
        disp('median')
        disp(exp(mu_ff1))
        disp(exp(mu_ff2))
        disp('std')
        disp(exp(sigma_ff1))
        disp(exp(sigma_ff2))
        disp(c);disp(n);
    end




    subplot(2,2,figOrder(regi))



    errorbar(binCenterff, meanCff, stdCff,'k');
    title(well_regList(regi))
    xlabel 'Length (ms)', ylabel 'Probability'
    xlim([100,10^(3.5)]);ylim([0,0.15]);
    yticks([0:0.02:0.14])
    xticks([10,100,1000,10000,100000])
    set(gca, 'xscale','log', 'fontsize',18)

    format_2x2_plot(regi)

    % Overlaying fit
    hold on
    fitvals = f_ff(log(binCenterff'));
    pf = plot(binCenterff, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';



end

saveas(gcf,strcat("Length_060623",'ff'),'png')

[~,~,stats1] = anova1(r',regLabel1ff);
stats1.means=stats1.means(~isnan(stats1.means));
stats1.n=stats1.n(stats1.n~=0);
[ff_c,ff_means]=multcompare(stats1,0.05,'on','','s');

xlabel ('log of mean')


%% Length FB

region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
regLabel1ff=[];regLabel1fb=[];
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
well_regList={"EC-DG","DG-CA3","CA3-CA1","CA1-EC"};
figOrder=[1 2 4 3];
if ~exist('graph', 'dir')
    mkdir('graph');
end

cd 'graph'
R_square=[];
r=[];

moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','mode','sd','n'});
figure( 'Position', [100 100 700 600])
binEdge = logspace(2,3.5,50); binCenter = convert_edges_2_centers(binEdge);
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];
    cross_anova_fb=[];
    c=0;
    n=0;
    for fi=1:9
        temp_setff=[];temp_setfb=[];
        tmp=tab.Len{1, fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff'
                temp_setff=[temp_setff,tmp{k+(regi-1)*5}.*0.04];
        
            elseif Index_fd{fi,k+(regi-1)*5}=='fb'
                temp_setfb=[temp_setfb,tmp{k+(regi-1)*5}.*0.04];
                if ~isempty(tmp{k+(regi-1)*5})
                    c=c+1;
                end
            end
        end
        totaln=numel(temp_setfb)+numel(temp_setff);
        [hist_objectfb,hist_plotfb,histcountfb(colifb,:)] = log_binned_histogram(temp_setfb, binEdge, "pdf",0,totaln);
        colifb=colifb+1;
        n=n+numel(temp_setfb);
        cross_anova_fb=[cross_anova_fb,temp_setfb];
        r=[r,temp_setfb];
        regLabel1fb = [regLabel1fb; repmat(categorical(region(regi)),numel(temp_setfb),1)];
    end
    
    histcountfb(~any(~isnan(histcountfb), 2),:)=[];
    meanCfb = mean(histcountfb,1); stdCfb = stdErr(histcountfb,1);
    meanCfb=meanCfb(~isnan(meanCfb));
    binCenterfb=binCenter(~isnan(meanCfb));
    Length_fb(regi)={cross_anova_fb};
    if regi==4
        [f_fb,gof_fb]=fit(log10(binCenterfb)', meanCfb', 'gauss1');%,'Lower',[0,0,0,0,2,0,0,3,0]);


        mu_fb(regi) = f_fb.b1; sigma_fb(regi) = f_fb.c1;
        [~,a2]=max(f_fb(log10(binCenterfb)));
        disp('mod')
        disp(binCenterfb(a2))
        t_val=f_fb(log10(binCenterfb))'.*binCenterfb;
        f_mean=mean(t_val);
        [~,a2]=min(abs(t_val-f_mean));
        disp('mean')
        disp(binCenterfb(a2))
        disp('median')
        disp(10^mu_fb(regi))
        disp("std")
        disp(10.^sigma_fb(regi))
        disp(c);disp(n);
    else
        [f_fb,gof_fb]=fit(log10(binCenterfb)', meanCfb', 'gauss2');
        mu_ff1 = f_fb.b1;
        mu_ff2 = f_fb.b2;
        sigma_ff1= f_fb.c1;
        sigma_ff2= f_fb.c2;
        disp('median')
        disp(10.^mu_ff1)
        disp(10.^mu_ff2)
        disp('std')
        disp(10.^sigma_ff1)
        disp(10.^sigma_ff2)
        disp(c);disp(n);
    end



    subplot(2,2,figOrder(regi))



    errorbar(binCenterfb, meanCfb, stdCfb,'k');
    title(well_regList(regi))
    xlabel 'Length (ms)', ylabel 'Probability'
    xlim([100,10^(3.5)]);ylim([0,0.14]);
    yticks([0:0.02:0.14])
    xticks([10,100,1000,10000,100000])
    set(gca, 'xscale','log', 'fontsize',18)

    format_2x2_plot(regi)

    % Overlaying fit
    hold on
    fitvals = f_fb(log10(binCenterfb)');
    pf = plot(binCenterfb, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    %    R_square(regi)=gof.rsquare;


end

saveas(gcf,strcat("Length_060623",'fb'),'png')

[~,~,stats1] = anova1(r',regLabel1fb);

[fb_c,fb_means]=multcompare(stats1,0.05,'on','','s');

cd ..
%% Length Paired ANOVA
for regi=1:4
if regi~=2
r=[Length_ff{regi},Length_ff{regi}];
n1=numel(Length_ff{regi});
n2=numel(Length_fb{regi});

regLabelPaired=[repmat(categorical(region(regi)),numel(n1),1);repmat(categorical(region(regi)),numel(n2),1)];

[~,~,stats1] = anova1(r',regLabelPaired);

[paired_c,paired_means]=multcompare(stats1,0.05,'on','','s');
end
end
%% Interval ff
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
regLabel1ff=[];regLabel1fb=[];
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
well_regList={"EC-DG","DG-CA3","CA3-CA1","CA1-EC"};
figOrder=[1 2 4 3];
if ~exist('graph', 'dir')
    mkdir('graph');
end

cd 'graph'
R_square=[];
r=[];
moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','mode','sd','n'});
figure( 'Position', [100 100 700 600])
binEdge = logspace(1,6,50); binCenter = convert_edges_2_centers(binEdge);
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];

    c=0;
    n=0;
    for fi=1:9
        temp_setff=[];temp_setfb=[];
        tmp=tab.In{1, fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff'
                temp_setff=[temp_setff,tmp{k+(regi-1)*5}.*0.04];
                if ~isempty(tmp{k+(regi-1)*5})
                    c=c+1;
                end
            elseif Index_fd{fi,k+(regi-1)*5}=='fb'
                temp_setfb=[temp_setfb,tmp{k+(regi-1)*5}.*0.04];

            end
        end
        totaln=numel(temp_setfb)+numel(temp_setff);
        [hist_objectff,hist_plotff,histcountff(coliff,:)] = log_binned_histogram(temp_setff, binEdge, "pdf",0,totaln);
        coliff=coliff+1;
        n=n+numel(temp_setff);
        r=[r,temp_setff];
        regLabel1ff = [regLabel1ff; repmat(categorical(region(regi)),numel(temp_setff),1)];
    end

    histcountff(~any(~isnan(histcountff), 2),:)=[];
    meanCff = mean(histcountff,1); stdCff = stdErr(histcountff,1);
    meanCff=meanCff(~isnan(meanCff));
    binCenterff=binCenter(~isnan(meanCff));


    if regi~=3 
        [f_ff,gof_ff]=fit(log(binCenterff'), meanCff', 'gauss1');
        mu_ff(regi) = f_ff.b1; sigma_ff(regi) = f_ff.c1/sqrt(2);%,'Lower',[0,0,0,0,2,0,0,3,0]);




        [~,a2]=max(f_ff(log(binCenterff)));
        disp('mod')
        disp(binCenterff(a2))
        t_val=f_ff(log(binCenterff))'.*binCenterff;
        f_mean=mean(t_val);
        [~,a2]=min(abs(t_val-f_mean));
        disp('mean')
        disp(binCenterff(a2))
        disp('median')
        disp(exp(mu_ff(regi)))
        disp('std')
        disp(exp(sigma_ff(regi)))

        disp(c);disp(n);

    else
        [f_ff,gof_ff]=fit(log(binCenterff'), meanCff', 'gauss2');
        mu_ff1 = f_ff.b1;
        mu_ff2 = f_ff.b2;
        sigma_ff1= f_ff.c1/sqrt(2);
        sigma_ff2= f_ff.c2/sqrt(2);
        disp('median')
        disp(exp(mu_ff1))
        disp(exp(mu_ff2))
        disp('std')
        disp(exp(sigma_ff1))
        disp(exp(sigma_ff2))
        disp(c);disp(n);
    end


    subplot(2,2,figOrder(regi))



    errorbar(binCenterff./1000, meanCff, stdCff,'k');
    title(well_regList(regi))
    xlabel 'Interval (s)', ylabel 'Probability'
    xlim([0.1,1000]);ylim([0,0.1]);
    yticks([0:0.01:0.1])
    xticks([1,10,100])
    set(gca, 'xscale','log', 'fontsize',18)

    format_2x2_plot(regi)

    % Overlaying fit
    hold on
    fitvals = f_ff(log(binCenterff'));
    pf = plot(binCenterff./1000, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';


end

saveas(gcf,strcat("interval_060623",'ff'),'png')
%
[~,~,stats1] = anova1(r',regLabel1ff);
stats1.means=stats1.means(~isnan(stats1.means));
stats1.n=stats1.n(stats1.n~=0);
[ff_c,ff_means]=multcompare(stats1,0.05,'on','','s');
% xlabel ('log of mean')
% figure
% % stats1.means(2)=-1;
% mu=exp(mu_ff(regi)+sigma_ff.^2/2);
%
% %bar(categorical(stats1.gnames'),mu)
% hold on;
%
% %  er=errorbar(categorical(stats1.gnames'),mu,sqrt(exp(sigma_ff.^2-1).*mu.^2));
% %   er.Color = [0 0 0];
% %  er.LineStyle = 'none';
%
% %ylim([0,4])
% ylabel ('mean(s)')
% saveas(gcf,strcat("interval anova",'ff'),'png')
cd ..
%% Interval fb

region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
regLabel1ff=[];regLabel1fb=[];
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
well_regList={"EC-DG","DG-CA3","CA3-CA1","CA1-EC"};
figOrder=[1 2 4 3];
if ~exist('graph', 'dir')
    mkdir('graph');
end

cd 'graph'
R_square=[];
r=[];

moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','mode','sd','n'});
figure( 'Position', [100 100 700 600])
binEdge = logspace(1,6,50); binCenter = convert_edges_2_centers(binEdge);
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];

    c=0;
    n=0;
    for fi=1:9
        temp_setff=[];temp_setfb=[];
        tmp=tab.In{1, fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff'
                temp_setff=[temp_setff,tmp{k+(regi-1)*5}.*0.04];

            elseif Index_fd{fi,k+(regi-1)*5}=='fb'
                temp_setfb=[temp_setfb,tmp{k+(regi-1)*5}.*0.04];
                if ~isempty(tmp{k+(regi-1)*5})
                    c=c+1;
                end
            end
        end
        totaln=numel(temp_setfb)+numel(temp_setff);
        [hist_objectfb,hist_plotfb,histcountfb(colifb,:)] = log_binned_histogram(temp_setfb, binEdge, "pdf",0,totaln);
        colifb=colifb+1;
        n=n+numel(temp_setfb);
        r=[r,temp_setfb];
        regLabel1fb = [regLabel1fb; repmat(categorical(region(regi)),numel(temp_setfb),1)];
    end

    histcountfb(~any(~isnan(histcountfb), 2),:)=[];
    meanCfb = mean(histcountfb,1); stdCfb = stdErr(histcountfb,1);
    meanCfb=meanCfb(~isnan(meanCfb));
    binCenterfb=binCenter(~isnan(meanCfb));
    


    if regi~=2
    [f_fb,gof_fb]=fit(log10(binCenterfb)', meanCfb', 'gauss1');%,'Lower',[0,0,0,0,2,0,0,3,0]);


    mu_fb(regi) = f_fb.b1; sigma_fb(regi) = f_fb.c1/sqrt(2);

    [~,a2]=max(f_fb(log10(binCenterfb)));
    disp('mod')
    disp(binCenterfb(a2))
    t_val=f_fb(log10(binCenterfb))'.*binCenterfb;
    f_mean=mean(t_val);
    [~,a2]=min(abs(t_val-f_mean));
    disp('mean')
    disp(binCenterfb(a2))
    disp('median')
    disp(10^mu_fb(regi))
    disp('std')
    disp(10.^sigma_fb(regi))
    disp(c);disp(n);
    else
        [f_fb,gof_fb]=fit(log10(binCenterfb)', meanCfb', 'gauss2');%,'Lower',[0,0,0,0,2,0,0,3,0]);


    mu_fb1 = f_fb.b1; sigma_fb(regi) = f_fb.c1;
    mu_fb2 = f_fb.b2;
    [~,a2]=max(f_fb(log10(binCenterfb)));
    disp('mod')
    disp(f_fb.b1);
    disp(f_fb.b2);
    disp(c);disp(n);
    end



    subplot(2,2,figOrder(regi))



    errorbar(binCenterfb./1000, meanCfb, stdCfb,'k');
    title(well_regList(regi))
    xlabel 'Interval (s)', ylabel 'Probability'
    xlim([0.1,1000]);ylim([0,0.1]);
    yticks([0:0.01:0.1])
    xticks([1,10,100])
    set(gca, 'xscale','log', 'fontsize',18)

    format_2x2_plot(regi)

    % Overlaying fit
    hold on
    fitvals = f_fb(log10(binCenterfb)');
    pf = plot(binCenterfb./1000, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    %    R_square(regi)=gof.rsquare;


end

saveas(gcf,strcat("interval_060623",'fb'),'png')

[~,~,stats2] = anova1(r',regLabel1fb);

[fb_c,fb_means]=multcompare(stats2,0.05,'on','','s');
% xlabel ('log of mean')
% figure
% %stats2.means(1)=-1;
% bar(categorical(stats2.gnames'),10.^stats2.means./1000)
% ylabel ('mean(s)')
% ylim([0,4])
% saveas(gcf,strcat("interval anova",'fb'),'png')
cd ..
%% feed direction comparison
% for regi=1:4
%     r=[r_ff{regi},r_fb{regi}];
%     regLabel = [repmat("FeedForward",numel(r_ff{regi}),1); repmat("Feedback",numel(r_fb{regi}),1)];
%     [~,~,stats] = anova1(log(r)',regLabel);
%     [fb_c,fb_means]=multcompare(stats,0.05,'on','','s');
% end





%% Amplitude ff
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
regLabel1ff=[];regLabel1fb=[];
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
well_regList={"EC-DG","DG-CA3","CA3-CA1","CA1-EC"};
figOrder=[1 2 4 3];
if ~exist('graph', 'dir')
    mkdir('graph');
end

cd 'graph'
R_square=[];
r=[];

moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','mode','sd','n'});
figure( 'Position', [100 100 700 600])
binEdge = logspace(0,3,50); binCenter = convert_edges_2_centers(binEdge);
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];

    c=0;
    n=0;
    for fi=1:9
        temp_setff=[];temp_setfb=[];
        tmp1=tab.IC{1, fi};
        tmp2=tab.Len{1,fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff'
                temp_setff=[temp_setff,tmp1{k+(regi-1)*5}./tmp2{k+(regi-1)*5}];
                if ~isempty(tmp1{k+(regi-1)*5})
                    c=c+1;
                end
            elseif Index_fd{fi,k+(regi-1)*5}=='fb'
                temp_setfb=[temp_setfb,tmp1{k+(regi-1)*5}./tmp2{k+(regi-1)*5}];
            end
        end
        totaln=numel(temp_setfb)+numel(temp_setff);
        [hist_objectff,hist_plotff,histcountff(coliff,:)] = log_binned_histogram(temp_setff, binEdge, "pdf",0,totaln);
        coliff=coliff+1;
        n=n+numel(temp_setff);
        r=[r,temp_setff];
        regLabel1ff = [regLabel1ff; repmat(categorical(region(regi)),numel(temp_setff),1)];
    end

    histcountff(~any(~isnan(histcountff), 2),:)=[];
    meanCff = mean(histcountff,1); stdCff = stdErr(histcountff,1);
    meanCff=meanCff(~isnan(meanCff));
    binCenterff=binCenter(~isnan(meanCff));



    if regi~=3  & regi~=2
        [f_ff,gof_ff]=fit(log10(binCenterff)', meanCff', 'gauss1');%,'Lower',[0,0,0,0,2,0,0,3,0]);


        mu_ff(regi) = f_ff.b1; sigma_ff(regi) = f_ff.c1/sqrt(2);

        [~,a2]=max(f_ff(log10(binCenterff)));
        disp('mod')
        disp(binCenterff(a2))
        t_val=f_ff(log10(binCenterff))'.*binCenterff;
        f_mean=mean(t_val);
        [~,a2]=min(abs(t_val-f_mean));
        disp('mean')
        disp(binCenterff(a2))
        disp('median')
        disp(10^mu_ff(regi))
        disp('std')
        disp(10.^sigma_ff(regi))
        disp(c);disp(n);
    elseif regi==3
        [f_ff,gof_ff]=fit(log10(binCenterff)', meanCff', 'gauss2');%,'Lower',[0,0,0,0,2,0,0,3,0]);


        mu_ff1 = f_ff.b1;  mu_ff2 = f_ff.b2; sigma_ff1 = f_ff.c1;sigma_ff2 = f_ff.c2;

        [~,a2]=max(f_ff(log10(binCenterff)));
        disp('mod')
        disp(binCenterff(a2))
        t_val=f_ff(log10(binCenterff))'.*binCenterff;
        f_mean=mean(t_val);
        [~,a2]=min(abs(t_val-f_mean));
        disp('mean')
        disp(binCenterff(a2))
        disp('median')
        disp(10^mu_ff1)
        disp(10^mu_ff2)
        disp('std')
        disp(10.^sigma_ff1)
        disp(10.^sigma_ff2)
        disp(c);disp(n);
    else
         [f_ff,gof_ff]=fit([1,2,3]', [0,0,0]', 'gauss1');%,'Lower',[0,0,0,0,2,0,0,3,0]);


        mu_ff(regi) = f_ff.b1; sigma_ff(regi) = f_ff.c1/sqrt(2);


    end


    subplot(2,2,figOrder(regi))


    if regi~=2
    errorbar(binCenterff, meanCff, stdCff,'k');
    title(well_regList(regi))
    xlabel 'Amplitude (uV)', ylabel 'Probability'
    xlim([0,10^(3)]);ylim([0,0.14]);
    yticks([0:0.02:0.14])
    xticks([0,10,100,1000])
    set(gca, 'xscale','log', 'fontsize',18)

    format_2x2_plot(regi)

    % Overlaying fit
    hold on
    fitvals = f_ff(log10(binCenterff)');
    pf = plot(binCenterff, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    %    R_square(regi)=gof.rsquare;
    end

end

saveas(gcf,strcat("amplitude_060623",'ff'),'png')

[~,~,stats1] = anova1(r',regLabel1ff);
stats1.means=stats1.means(~isnan(stats1.means));
stats1.n=stats1.n(stats1.n~=0);
[ff_c,ff_means]=multcompare(stats1,0.05,'on','','s');
% xlabel ('log of mean')
% figure
% %stats1.means(2)=-1;
% bar(categorical(stats1.gnames'),10.^stats1.means)
% ylabel ('mean(uV)')
% ylim([0, 45])
% saveas(gcf,strcat("Amplitude anova",'ff'),'png')
cd ..

%% Amplitude fb

region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
regLabel1ff=[];regLabel1fb=[];
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
well_regList={"EC-DG","DG-CA3","CA3-CA1","CA1-EC"};
figOrder=[1 2 4 3];
if ~exist('graph', 'dir')
    mkdir('graph');
end

cd 'graph'
R_square=[];
r=[];

moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','mode','sd','n'});
figure( 'Position', [100 100 700 600])
binEdge = logspace(0,3,50); binCenter = convert_edges_2_centers(binEdge);
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];

    c=0;
    n=0;
    for fi=1:9
        temp_setff=[];temp_setfb=[];
        tmp1=tab.IC{1, fi};
        tmp2=tab.Len{1,fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff'
                temp_setff=[temp_setfb,tmp1{k+(regi-1)*5}./tmp2{k+(regi-1)*5}];

            elseif Index_fd{fi,k+(regi-1)*5}=='fb'
                temp_setfb=[temp_setfb,tmp1{k+(regi-1)*5}./tmp2{k+(regi-1)*5}];
                if ~isempty(tmp2{k+(regi-1)*5})
                    c=c+1;
                end
            end
        end
        totaln=numel(temp_setfb)+numel(temp_setff);

        [hist_objectfb,hist_plotfb,histcountfb(colifb,:)] = log_binned_histogram(temp_setfb, binEdge, "pdf",0,totaln);
        colifb=colifb+1;
        n=n+numel(temp_setfb);
        r=[r,temp_setfb];
        regLabel1fb = [regLabel1fb; repmat(categorical(region(regi)),numel(temp_setfb),1)];
    end

    histcountfb(~any(~isnan(histcountfb), 2),:)=[];
    meanCfb = mean(histcountfb,1); stdCfb = stdErr(histcountfb,1);
    meanCfb=meanCfb(~isnan(meanCfb));
    binCenterfb=binCenter(~isnan(meanCfb));

    if regi~=1 & regi~=2
        [f_fb,gof_fb]=fit(log10(binCenterfb)', meanCfb', 'gauss1');%,'Lower',[0,0,0,0,2,0,0,3,0]);
        mu_fb(regi) = f_fb.b1; sigma_fb(regi) = f_fb.c1;
        disp('median')
        disp(10^mu_fb(regi))
        disp('std')
        disp(10.^sigma_fb(regi))
        disp(c);disp(n);
    elseif regi~=2
        [f_fb,gof_fb]=fit(log10(binCenterfb)', meanCfb', 'gauss2');%,'Lower',[0,0,0,0,2,0,0,3,0]);
        mu_fb(regi) = f_fb.b1; sigma_fb1 = f_fb.c1; sigma_fb2 = f_fb.c2;
        disp('median')
        disp(exp(f_fb.b1))
        disp(exp(f_fb.b2))
        disp('std')
        disp(exp(sigma_fb1))
        disp(exp(sigma_fb2))
        disp(c);disp(n);
    elseif regi==2
        [f_fb,gof_fb]=fit(log10(binCenterfb)', meanCfb', 'gauss3','Lower',[0,1,0,0,0,0,0,2,0.08],'upper',[Inf,1.36,Inf,Inf,Inf,Inf,Inf,3,Inf]);
        mu_fb(regi) = f_fb.b1; sigma_fb1 = f_fb.c1; sigma_fb2 = f_fb.c2;
        disp('median')
        disp(10^(f_fb.b1))
        disp(10^(f_fb.b2))
         disp(10^(f_fb.b3))
        disp(c);disp(n);
    end



    subplot(2,2,figOrder(regi))



    errorbar(binCenterfb, meanCfb, stdCfb,'k');
    title(well_regList(regi))
    xlabel 'Amplitude (uV)', ylabel 'Probability'
    xlim([0,10^(3)]);ylim([0,0.14]);
    yticks([0:0.02:0.14])
    xticks([0,10,100,1000])
    set(gca, 'xscale','log', 'fontsize',18)

    format_2x2_plot(regi)

    % Overlaying fit
    hold on
    fitvals = f_fb(log10(binCenterfb)');
    pf = plot(binCenterfb, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    %    R_square(regi)=gof.rsquare;


end

saveas(gcf,strcat("amplitude_060623",'fb'),'png')
[~,~,stats1] = anova1(r',regLabel1fb);

[fb_c,fb_means]=multcompare(stats1,0.05,'on','','s');
% xlabel ('log of mean')
% stats1.means(1)=-1;
%
% figure
% bar(categorical(stats1.gnames'),10.^stats1.means)
% ylabel ('mean(uV)')
% ylim([0, 45])
% saveas(gcf,strcat("Amplitude anova",'fb'),'png')
cd ..




%% feed direction comparison
% for regi=1:4
%     r=[r_ff{regi},r_fb{regi}];
%     regLabel = [repmat("FeedForward",numel(r_ff{regi}),1); repmat("Feedback",numel(r_fb{regi}),1)];
%     [~,~,stats1] = anova1(log(r)',regLabel);
%     [fb_c,fb_means]=multcompare(stats1,0.05,'on','','s');
% end
%
%
%
%


















%% functions

function [ hist_object,hist_plot, hist_count] = log_binned_histogram( vec, binning_prarmeter, norm_method, if_plot,totaln )
%Populates positive log binned histogram and return histogram handle and fig
%binning parameter either no. of bins or Bin edges

if ~exist('norm_method','var'), norm_method = "none"; end
if ~exist('if_plot','var'), if_plot=1; end

hist_plot=[];
vec(vec <= 0)=[];
if length(binning_prarmeter)==1
    xEdge = logspace(round(log10(min(vec))), round(log10(max(vec))), binning_prarmeter);
else
    xEdge = binning_prarmeter;
end
xCenter = convert_edges_2_centers(xEdge);
hist_object=  histogram(vec, xEdge,'Visible','off');
if norm_method == "pdf"
    hist_count = hist_object.BinCounts/ totaln;
elseif norm_method == "cdf"
    hist_count = fliplr(cumsum(fliplr(hist_object.BinCounts)))./sum(hist_object.BinCounts);
elseif norm_method == "count-cdf"
    hist_count = fliplr(cumsum(fliplr(hist_object.BinCounts)));
else
    hist_count = hist_object.BinCounts;
end

if(if_plot)
    hist_plot=loglog( xCenter, hist_count,'LineWidth',2);
    set(gca,'XScale','log')

end

end
function [ centers ] = convert_edges_2_centers( edges  )
%Converts bin-edges to bin-centers
centers = conv(edges, [0.5 0.5], 'valid');
end
function format_2x2_plot( i)

switch i
    case 1
        xlabel '';
    case 2
        xlabel ''; ylabel '';
    case 3
        ylabel '';
end

end