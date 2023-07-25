%% spiking fitting
close all
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
well_regList={"EC-DG","DG-CA3","CA3-CA1","CA1-EC"};
figOrder=[1 2 4 3];
if ~exist('graph', 'dir')
    mkdir('graph');
end
s=[];
cd 'graph'
for i=1:4
    figure
    co=1;
    for j=1:6
        tmp=tab.pc{1, j};
        temp_set=[tmp{1+(i-1)*5},tmp{2+(i-1)*5},tmp{3+(i-1)*5},tmp{4+(i-1)*5},tmp{5+(i-1)*5}];
        s=[s,temp_set];
    end

    histogram(s);
    title(strcat("%Spike Positive",region(i)));
    xlabel("% positive")
    saveas(gcf,strcat("Spike Posiiton",region(i)),'png')
end

%% Spiking position
for i=1:4
    figure
    co=1;
    for j=1:6
        tmp1=tab.pl{1, j};
        tmp2=tab.ph{1, j};
        tmp3=tab.data3{1, j};
        tmp4=[];
        for k=1:5
            temp_set1=tmp1{k+(i-1)*5}./25;
            temp_set2=tmp2{k+(i-1)*5}./25;
            temp_set3=tmp3{k+(i-1)*5};
            for l=1:numel(temp_set1)
                con(l)=numel(temp_set3(temp_set3<=temp_set2(l) & temp_set3>=temp_set1(l)));
            end
            tmp4=[tmp4,con];
        end
        binEdge = logspace(-1,4,100); binCenter = convert_edges_2_centers(binEdge);
        [ hist_object, hist_plot, prob(co,:)] =log_binned_histogram(tmp4,binEdge,'pdf',0);
        co=co+1;
    end
    prob(~any(~isnan(prob),2), :)=[];
    meanC=mean(prob,1); stdC=stdErr(prob,1);
    [f,gof]=fit(log(binCenter)', meanC', 'gauss1');
    plot(f,log(binCenter), meanC)
    R_square(i)=gof.rsquare;
    %legend("EC-DG","","DG-CA3","","CA3-CA1","","CA1-EC");
    xlabel("ln(mumber)")
    ylim([0 0.05]);
    title(strcat("Spike Count within Spindle ",region(i)));
    ylabel("Probability")
    %xlim([1e3,1e5])
    hold off;
    saveas(gcf,strcat("Spike Count",region(i)),'png')
end


moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','mode','sd','n'});
figure( 'Position', [100 100 700 600])
%figure( 'Position', [100 100 2800 600])
binEdge = logspace(-1,4,100); binCenter = convert_edges_2_centers(binEdge);
co=0;
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:6
        tmp1=tab.pl{1, fi};
        tmp2=tab.ph{1, fi};
        tmp3=tab.data3{1, fi};
        tmp4=[];
        for k=1:5
            temp_set1=tmp1{k+(regi-1)*5}./25;
            temp_set2=tmp2{k+(regi-1)*5}./25;
            temp_set3=tmp3{k+(regi-1)*5};
            for l=1:numel(temp_set1)
                con(l)=numel(temp_set3(temp_set3<=temp_set2(l) & temp_set3>=temp_set1(l)));
            end
            tmp4=[tmp4,con];
        end
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(tmp4, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);

    % calculating mean and sd

    [f,gof]=fit(log(binCenter)', meanC', 'gauss1');

    %     if regi==3
    %         f = fit(log(binCenter)', meanC', 'gauss1','Exclude', find(ismember(meanC',maxk(meanC',2))));
    %     end

    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi+co)=regi;

    %fitted
    moment_table.mean(regi+co) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi+co) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi+co) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi+co) ^2));
    moment_table.n(regi+co)=length(f);

    % plotting data
    errorbar(binCenter, meanC, stdC,'k');
    title(well_regList(regi))
    xlabel 'number of spike', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18)
    xlim([1e-1 10^(4)]), ylim([0 0.05])
    format_2x2_plot(regi)

    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';

    %     moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    %     moment_table.emp_sd(regi) = sqrt(var(binCenter,fitvals));

end
saveas(gcf,"Spike count real",'png')
moment_table
cd ..

%% Correlation
if ~exist('graph', 'dir')
    mkdir('graph');
end
cw=load('matching_table_cw.mat');
ccw=load('matching_table_ccw.mat');
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
cd 'graph'

for i=1:4
    figure
    co=1;
    temp_set1=[]; temp_set2=[];
    for j=1:6
        tmp1=tab.rspin{1, j};
        temp_set1=[temp_set1,tmp1{1+(i-1)*5},tmp1{2+(i-1)*5},tmp1{3+(i-1)*5},tmp1{4+(i-1)*5},tmp1{5+(i-1)*5}];

        tmp2=tab.rspike{1, j};
        temp_set2=[temp_set2,tmp2{1+(i-1)*5},tmp2{2+(i-1)*5},tmp2{3+(i-1)*5},tmp2{4+(i-1)*5},tmp2{5+(i-1)*5}];


    end

    [cor,lags]=xcorr(temp_set1,temp_set2,400,'normalized');
    Peak=max(cor);
    PeakIdx=find(cor==Peak);
    PP=lags(PeakIdx);
    stem(lags,cor)
    text(lags(PeakIdx), Peak+0.1, strcat(num2str(PP)," , ",num2str(Peak)))
    ylabel("Correlation of Spike train With Spindle")
    title(region(i))
    ylim([0 1])
    xline(0)
    xlabel("Spindle before Spike (ms) Spike before Spindle")
    saveas(gcf,strcat("Corr",region(i)),'png')

end
if ~exist('xcorr', 'dir')
    mkdir('xcorr');
end
cd 'xcorr'
qcw=[];
for i=1:9
    figure

    co=1;
    tmp1=tab.rspin{1, i};
    tmp2=tab.rspike{1, i};
    if ~exist(s_folders(i), 'dir')
        mkdir(s_folders(i));
    end
    cd(s_folders(i))
    if is_cw(i)==1
        na=table2cell(cw.matching_table  );
        qcw=[qcw,'T'];
    else
        na=table2cell(ccw.matching_table  );
        qcw=[qcw,'F'];
    end
    for j=1:20
        temp_set1=tmp1{j};
        temp_set2=tmp2{j};
        if ~isempty(temp_set2)

            [cor,lags]=xcorr(temp_set1,temp_set2,400,'normalized');
            Peak=max(cor);
            PeakIdx=find(cor==Peak);
            stem(lags,cor)
            PP=lags(PeakIdx);
            text(lags(PeakIdx), Peak+0.1, strcat(num2str(PP)," , ",num2str(Peak)))
            title(na{j,2})
            xline(0)
            xlabel("Spindle before Spike (ms) Spike before spindle")
            ylabel("Correlation of Spike train With Spindle")
            ylim([0 1])
            saveas(gcf,strcat("Corr",na{j,2}),'png')
        end
    end
    cd ..
    close all

end


cd ../..

%% polar histogram for phase angle

cutoff=[25,15,15,15];
set(0,'defaultAxesFontSize',20)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
slope=[];
Name_test=[];
if ~exist('graph', 'dir')
    mkdir('graph');
end
cw=load('matching_table_cw.mat');
ccw=load('matching_table_ccw.mat');
region=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
cd 'graph'
for i=1:4

    co=1;
    temp_set1=[];
    for j=1:9
        tmp1=tab.angle{1, j};
        temp_set1=[temp_set1,tmp1{1+(i-1)*5},tmp1{2+(i-1)*5},tmp1{3+(i-1)*5},tmp1{4+(i-1)*5},tmp1{5+(i-1)*5}];
    end
%     polarhistogram(temp_set1,50)
%    [avg_spike_angle(i),magnitude(i),~,~]=circ_mean(temp_set1);
    % hold on;
    % polarplot(avg_spike_angle(i),magnitude(i),'r-')
%     title (region(i))
%     saveas(gcf,strcat("Angle",region(i)),'png')

end
if ~exist('angle', 'dir')
    mkdir('angle');
end
cd 'angle'
for i=1:4
    co=1;
    temp_set1=[];
    reg_setff_h=[];
    reg_setff_l=[];
    reg_setfb_h=[];
    reg_setfb_l=[];
    for j=1:9
        if ~exist(s_folders(j), 'dir')
            mkdir(s_folders(j));
        end
        cd(s_folders(j))
        tmp1=tab.angle{1, j};
        tmp2=tab.IC{1,j};
        tmp3=tab.Len{1,j};
        tmp4=tab.spik_n{1,j};
        if is_cw(j)==1
            tb=table2cell(cw.matching_table);
        else
            tb=table2cell(ccw.matching_table);
        end
        for k=1:5
            amp_a=tmp2{k+(i-1)*5}./tmp3{k+(i-1)*5};
            tmp_ind=[];
            for l=1:numel(tmp4{k+(i-1)*5})
               % if amp_a(l)>=cutoff(i)
                    index_g=ones(1,tmp4{k+(i-1)*5}(l));
                %else
                  %  index_g=zeros(1,tmp4{k+(i-1)*5}(l));
                %end
                tmp_ind=[tmp_ind,index_g];
            end
            angle_plot=tmp1{k+(i-1)*5};





            if Index_fd{j, k+(i-1)*5} == 'ff'

%                figure
%                 h=polarhistogram(angle_plot(tmp_ind==1),50,'FaceColor','r');
% 
%                 title (strcat(region(i)," - ",tb{k+(i-1)*5,2}));
%                 hold on;
%                [avg_spike_angle,magnitude,~,~]=circ_mean(angle_plot(tmp_ind==1));
%                 polarplot([0, avg_spike_angle], [0,magnitude*max(h.Values)], '-r', 'LineWidth',3)
%                 hold off;
%                 saveas(gcf,strcat(region(i)," - ",tb{k+(i-1)*5,2}),'png')
                reg_setff_h=[reg_setff_h,angle_plot(tmp_ind==1)];

%                 figure
%                 h=polarhistogram(angle_plot(tmp_ind==0),50,'FaceColor','m');
%                 title (strcat(region(i)," - ",tb{k+(i-1)*5,2}));
%                 hold on;
%                 [avg_spike_angle,magnitude,~,~]=circ_mean(angle_plot(tmp_ind==0));
%                 polarplot([0, avg_spike_angle], [0,magnitude*max(h.Values)], '-r', 'LineWidth',3)
%                 hold off;
%                 saveas(gcf,strcat(region(i)," - ",tb{k+(i-1)*5,2}," smaller"),'png')
%                 reg_setff_l=[reg_setff_l,angle_plot(tmp_ind==0)];

            end

            if Index_fd{j, k+(i-1)*5} == 'fb'
%                figure
%                 h=polarhistogram(angle_plot(tmp_ind==1),50,'FaceColor','b');
%                 title (strcat(region(i)," - ",tb{k+(i-1)*5,2}));
%                 hold on;
%               %  [avg_spike_angle,magnitude,~,~]=circ_mean(angle_plot(tmp_ind==1));
%                 polarplot([0, avg_spike_angle], [0,magnitude*max(h.Values)], '-r', 'LineWidth',3)
%                 hold off;
%                 saveas(gcf,strcat(region(i)," - ",tb{k+(i-1)*5,2}),'png')
                reg_setfb_h=[reg_setfb_h,angle_plot(tmp_ind==1)];



%                 figure
%                 h=polarhistogram(angle_plot(tmp_ind==0),50,'FaceColor','c');
%                 title (strcat(region(i)," - ",tb{k+(i-1)*5,2}));
%                 hold on;
%                 [avg_spike_angle,magnitude,~,~]=circ_mean(angle_plot(tmp_ind==0));
%                 polarplot([0, avg_spike_angle], [0,magnitude*max(h.Values)], '-r', 'LineWidth',3)
%                 hold off;
%                 saveas(gcf,strcat(region(i)," - ",tb{k+(i-1)*5,2}," smaller"),'png')
%                 reg_setfb_l=[reg_setfb_l,angle_plot(tmp_ind==0)];
            end
        end
        cd ..







    end

    %[avg_spike_angle(i),magnitude(i),~,~]=circ_mean(temp_set1);
    % hold on;
    % polarplot(avg_spike_angle(i),magnitude(i),'r-')

    % FF h

    sz=0.0175;
figure
%     if i~=3
%     figure
%         ax = polaraxes;
    h=polarhistogram( reg_setff_h,50,'FaceColor','r');
    hold on;
    bc=[];bc1=[];bc2=[];
if i~=2 & i~=3
bc=convert_edges_2_centers(h.BinEdges);
        %bc(bc<=0)=2*pi+bc(bc<=0);
            [~,I]=min(h.BinCounts);
    bc1=[bc(I:end),bc(1:I-1)+(2.*pi)];
       bc2=[h.BinCounts(I:end),h.BinCounts(1:I-1)];
%     [pval_ff(i), m] = circ_otest(bc, sz,h.BinCounts);
%     disp(pval_ff(i))
%     disp(m)
hold off
disp('ff')
disp(region(i))
    f=fit(bc1',bc2'./sum(bc2'),'poly1');
    [~,I1]=max(bc2);[~,I2]=min(bc2);
    effect_size=(bc1(I1)-bc1(I2))/std(bc2);
    slope=[slope,f.p1];
            stats = regstats(bc1',bc2'./sum(bc2'),'linear');
        disp(effect_size);

figure
plot(bc1,bc2'./sum(bc2'),'o')
hold on
plot(bc1,f(bc1))
title (strcat("FF ",region(i)))
Name_test=[Name_test,strcat(region(i)," ff")];
ylim([0.01,0.03])
xlabel('angle(radians)')
ylabel("% of total")
hold off
elseif i==3
bc=convert_edges_2_centers(h.BinEdges);
        %bc(bc<=0)=2*pi+bc(bc<=0);
            [~,I]=min(h.BinCounts);
    bc1=[bc(I:end),bc(1:I-1)+(2.*pi)];
       bc2=[h.BinCounts(I:end),h.BinCounts(1:I-1)];
       bc1_1=bc1(bc1<=1.19);bc1_2=bc1(bc1>1.19 );
       bc2_1=bc2(bc1<=1.19);bc2_2=bc2(bc1>1.19 );
%     [pval_ff(i), m] = circ_otest(bc, sz,h.BinCounts);
%     disp(pval_ff(i))
%     disp(m)
hold off
disp('ff')
disp(region(i))
    f1=fit(bc1_1',bc2_1'./sum(bc2'),'poly1');
    f2=fit(bc1_2',bc2_2'./sum(bc2'),'poly1');
    slope=[slope,f1.p1.f2.p1];
figure
plot(bc1,bc2'./sum(bc2'),'o')
hold on
plot(bc1_1,f1(bc1_1))
plot(bc1_2,f2(bc1_2))
title (strcat("FF ",region(i)))
Name_test=[Name_test,strcat(region(i)," ff")];
ylim([0.01,0.03])
xlabel('angle(radius)')
ylabel("% of total")
hold off
        stats1 = regstats(bc1_1',bc2_1'./sum(bc2'),'linear');
            [~,I1]=max(bc2_1);[~,I2]=min(bc2_1);
    effect_size=(bc1(I1)-bc1(I2))/std(bc2_1);


        disp(effect_size);
        stats1 = regstats(bc1_2',bc2_2'./sum(bc2'),'linear');
            [~,I1]=max(bc2_2);[~,I2]=min(bc2_2);
    effect_size=(bc1(I1)-bc1(I2))/std(bc2_2);

        disp(effect_size);
else
    slope=[slope,NaN];
    Name_test=[Name_test,strcat(region(i)," ff")];
end
%     ax.FontSize = 20;
%     title (strcat(region(i)," ff "));
%     hold on;
%     [avg_spike_angle,magnitude,~,~]=circ_mean(reg_setff_h);

%     polarplot([0, avg_spike_angle], [0,magnitude*sum(h.Values)/50], '-g', 'LineWidth',3)

    hold off;
    saveas(gcf,strcat(region(i)," ff "),'png')
% 
%     else
%     figure
%     ax = polaraxes;
%     h=polarhistogram( reg_setff_h,50,'FaceColor','r');
%     ax.FontSize = 20;
%     title (strcat(region(i)," ff "));
%     hold on;
% %     [avg_spike_angle1,magnitude1,~,~]=circ_mean(reg_setff_h(reg_setff_h<=0.524 & reg_setff_h>=-1.047));
% % [avg_spike_angle2,magnitude2,~,~]=circ_mean(reg_setff_h(reg_setff_h<=3.141 & reg_setff_h>=1.471));
% 
%     polarplot([0, avg_spike_angle1], [0,magnitude1*numel(reg_setff_h(reg_setff_h<=0.524 & reg_setff_h>=-1.047))/12], '-g', 'LineWidth',3)
% polarplot([0, avg_spike_angle2], [0,magnitude2*numel(reg_setff_h(reg_setff_h<=3.141 & reg_setff_h>=1.471))/12], '-g', 'LineWidth',3)
%     hold off;
%     saveas(gcf,strcat(region(i)," ff "),'png')
%     bc=convert_edges_2_centers(h.BinEdges);
%         bc(bc<=0)=2*pi+bc(bc<=0);
%     [pval_ff(i), m] = circ_otest(bc,sz,h.BinCounts);
%     disp(pval_ff(i))
%     disp(m)
%          f=fit([bc,h.BinCounts],'poly1')
%     end

    % FF s
%     figure
%     h=polarhistogram( reg_setff_l,50,'FaceColor','m');
% 
%     title (strcat(region(i)," ff smaller"));
%     hold on;
%     [avg_spike_angle,magnitude,~,~]=circ_mean( reg_setff_l);
%     polarplot([0, avg_spike_angle], [0,magnitude*max(h.Values)], '-r', 'LineWidth',3)
%     hold off;
%     saveas(gcf,strcat(region(i)," ff smaller"),'png')
    % FB h
%     figure
figure
    ax = polaraxes;
    h=polarhistogram( reg_setfb_h,50,'FaceColor','b');
    hold on;
    ax.FontSize = 20;
    title (strcat(region(i)," fb"));
    Name_test=[Name_test,strcat(region(i)," fb")];
    bc=[];bc1=[];bc2=[];
%     hold on;
%     [avg_spike_angle,magnitude,~,~]=circ_mean(reg_setfb_h);
%    polarplot([0, avg_spike_angle], [0,magnitude*sum(h.Values)/50], '-g', 'LineWidth',3)
%    hold off;
%    saveas(gcf,strcat(region(i)," fb"),'png')
    bc=convert_edges_2_centers(h.BinEdges);
  %  bc(bc<=0)=2*pi+bc(bc<=0);
    [~,I]=min(h.BinCounts);
    bc1=[bc(I:end),bc(1:I-1)+(2.*pi)];
    bc2=[h.BinCounts(I:end),h.BinCounts(1:I-1)];
    hold off;
%     [pval_fb(i), m] = circ_otest(bc,sz,h.BinCounts);
%     disp(pval_fb(i))
%     disp(m)
disp('fb')
disp(region(i))
    f=fit(bc1',bc2'./sum(bc2'),'poly1');
        slope=[slope,f.p1];
            [~,I1]=max(bc2);[~,I2]=min(bc2);
    effect_size=(bc1(I1)-bc1(I2))/std(bc2);
        stats = regstats(bc1',bc2'./sum(bc2'),'linear');
        disp(effect_size);
    figure
plot(bc1,bc2./sum(bc2'),'o')
hold on
ylim([0.01,0.03])
plot(bc1,f(bc1))
title (strcat("FB ",region(i)))
xlabel('angle(radians)')
ylabel("% of total")
hold off
    saveas(gcf,strcat(region(i)," fb "),'png')
    % FB s
%     figure
%     h=polarhistogram( reg_setfb_l,50,'FaceColor','c');
% 
%     title (strcat(region(i)," fb smaller"));
%     hold on;
%     [avg_spike_angle,magnitude,~,~]=circ_mean(reg_setfb_l);
%     polarplot([0, avg_spike_angle], [0,magnitude*max(h.Values)], '-r', 'LineWidth',3)
%     hold off;
%     saveas(gcf,strcat(region(i)," fb smaller"),'png')
end
% 
% Names = categorical(Name_test);
%   [~,~,stats] = anova1(slope);
% [~,~]=multcompare(stats,0.05,'on','','s');


cd ..

cd ..
%% FUNCTION

function [ hist_object,hist_plot, hist_count] = log_binned_histogram( vec, binning_prarmeter, norm_method, if_plot )
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
    hist_count = hist_object.BinCounts/ sum(hist_object.BinCounts);
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