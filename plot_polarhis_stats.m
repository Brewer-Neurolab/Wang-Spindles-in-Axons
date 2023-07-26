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

                    index_g=ones(1,tmp4{k+(i-1)*5}(l));


                tmp_ind=[tmp_ind,index_g];
            end
            angle_plot=tmp1{k+(i-1)*5};





            if Index_fd{j, k+(i-1)*5} == 'ff'

                reg_setff_h=[reg_setff_h,angle_plot(tmp_ind==1)];

            end

            if Index_fd{j, k+(i-1)*5} == 'fb'

                reg_setfb_h=[reg_setfb_h,angle_plot(tmp_ind==1)];

            end
        end
        cd ..

    end



    sz=0.0175;
figure

    h=polarhistogram( reg_setff_h,50,'FaceColor','r');
    hold on;
    bc=[];bc1=[];bc2=[];
if i~=2 & i~=3
bc=convert_edges_2_centers(h.BinEdges);

            [~,I]=min(h.BinCounts);
    bc1=[bc(I:end),bc(1:I-1)+(2.*pi)];
       bc2=[h.BinCounts(I:end),h.BinCounts(1:I-1)];

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

            [~,I]=min(h.BinCounts);
    bc1=[bc(I:end),bc(1:I-1)+(2.*pi)];
       bc2=[h.BinCounts(I:end),h.BinCounts(1:I-1)];
       bc1_1=bc1(bc1<=1.19);bc1_2=bc1(bc1>1.19 );
       bc2_1=bc2(bc1<=1.19);bc2_2=bc2(bc1>1.19 );

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
    hold off;
    saveas(gcf,strcat(region(i)," ff "),'png')

figure
    ax = polaraxes;
    h=polarhistogram( reg_setfb_h,50,'FaceColor','b');
    hold on;
    ax.FontSize = 20;
    title (strcat(region(i)," fb"));
    Name_test=[Name_test,strcat(region(i)," fb")];
    bc=[];bc1=[];bc2=[];

    bc=convert_edges_2_centers(h.BinEdges);

    [~,I]=min(h.BinCounts);
    bc1=[bc(I:end),bc(1:I-1)+(2.*pi)];
    bc2=[h.BinCounts(I:end),h.BinCounts(1:I-1)];
    hold off;

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
  
end

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