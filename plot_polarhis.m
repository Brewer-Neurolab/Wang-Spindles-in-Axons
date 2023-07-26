%% polar histogram for phase angle

cutoff=[25,15,15,15];
set(0,'defaultAxesFontSize',20)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)

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


    sz=0.0175; % Bin size for polar histogram
    % FF histogram
    if i~=3
        figure
        ax = polaraxes;
        h=polarhistogram( reg_setff_h,50,'FaceColor','r');
        bc=convert_edges_2_centers(h.BinEdges);
        bc(bc<=0)=2*pi+bc(bc<=0);
        [pval_ff(i), m] = circ_otest(bc, sz,h.BinCounts);
        disp(pval_ff(i))
        disp(m)

        ax.FontSize = 20;
        title (strcat(region(i)," ff "));
        hold on;
        [avg_spike_angle,magnitude,~,~]=circ_mean(bc,h.BinCounts);

        polarplot([0, avg_spike_angle], [0,magnitude], '-g', 'LineWidth',3) %*sum(h.Values)/50

        hold off;
        saveas(gcf,strcat(region(i)," ff "),'png')
        disp(sum(h.BinCounts))
    else
        figure
        ax = polaraxes;
        h=polarhistogram( reg_setff_h,50,'FaceColor','r');
        [ centers ] = convert_edges_2_centers( h.BinEdges  );
        ax.FontSize = 20;
        title (strcat(region(i)," ff "));
        hold on;
        [avg_spike_angle1,magnitude1,~,~]=circ_mean(centers(centers<=0.524 & centers>=-1.047),h.BinCounts((centers<=0.524 & centers>=-1.047)));
        [avg_spike_angle2,magnitude2,~,~]=circ_mean(centers(centers<=3.141 & centers>=1.471),h.BinCounts((centers<=3.141 & centers>=1.471)));
        [avg_spike_angle3,magnitude3,~,~]=circ_mean(centers,h.BinCounts);
        polarplot([0, avg_spike_angle1], [0,magnitude1], '-y', 'LineWidth',3)
        polarplot([0, avg_spike_angle2], [0,magnitude2], '-y', 'LineWidth',3)
        polarplot([0, avg_spike_angle3], [0,magnitude3], '-g', 'LineWidth',3)
        hold off;
        saveas(gcf,strcat(region(i)," ff "),'png')
        bc=convert_edges_2_centers(h.BinEdges);
        bc(bc<=0)=2*pi+bc(bc<=0);
        [pval_ff(i), m] = circ_otest(bc,sz,h.BinCounts);
        disp(pval_ff(i))
        disp(m)
        disp(sum(h.BinCounts))
    end


    % FB histogram
    figure
    ax = polaraxes;
    h=polarhistogram( reg_setfb_h,50,'FaceColor','b');
    [ centers ] = convert_edges_2_centers( h.BinEdges  );
    ax.FontSize = 20;
    title (strcat(region(i)," fb"));
    hold on;
    [avg_spike_angle,magnitude,~,~]=circ_mean(centers,h.BinCounts);
    polarplot([0, avg_spike_angle], [0,magnitude], '-g', 'LineWidth',3)
    hold off;
    saveas(gcf,strcat(region(i)," fb"),'png')
    bc=convert_edges_2_centers(h.BinEdges);
    bc(bc<=0)=2*pi+bc(bc<=0);
    [pval_fb(i), m] = circ_otest(bc,sz,h.BinCounts);
    disp(pval_fb(i))
    disp(m)
    disp(sum(h.BinCounts))

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