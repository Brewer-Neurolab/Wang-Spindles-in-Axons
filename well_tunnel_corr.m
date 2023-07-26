%% data Cal for DG-CA3
set(0,'defaultAxesFontSize',32)
dir_lis=["H8","H9","H10","H11","H12","J8","J9","J10","J11","J12","K8","K9","K10","K11","L8","L9","L10","M8","M9"];

load('test 32 CA3 400ms thr.mat')
wells=tab.rspin;


load('test 31 400ms thr.mat')
tunnel=tab.rspin;


is_cw=[1,1,0,1,0,1,0,1,0];
load('Index_fd_CA3_V2.mat');
for fi=1:9
    wells_temp=wells{1,fi};
    tunnel_temp=tunnel{1,fi};
    for i=6:15
        tunnel_corr(i-5)={tunnel_temp{1,i}};
    end
    if is_cw(fi)==1
        load matching_table_cw.mat
    else
        load matching_table_ccw.mat
    end
    c=0;
    Bestcorr=[];
    Bestlagg=[];
    fd_dr={};
    fd_dr_best={};
    all_lag={};
    all_cor={};
    all_l_std={};
    all_c_std={};
    all_n={};
    all_name={};
    best_name={};
    color=[];

    for j=1:19
        if ~isempty(wells_temp{1,j})
            
            wells_corr=wells_temp{1,j};
            corr_r=[];
            temp_lag=[];
            names={};
            temp_fd={};
            temp_lag_std=[];
            temp_cor_std=[];
            temp_n=[];
            c_k=0;

            for k=1:5 %1:5
                if ~isempty(tunnel_corr{k})
                    c_k=c_k+1;
                    temp_r=[];
                    temp_l=[];
                    c_m=0;
                    %                     for m=1:numel(tab.ph{1, fi}{1, k+5})
                    %                         temp_low=floor(tab.pl{1, fi}{1, k+5}./25);
                    %                         temp_high=floor(tab.ph{1, fi}{1, k+5}./25);
                    %                         temp_low(temp_low==0)=1;
                    %                         if temp_high(m)-temp_low(m)>=10
                    %                         %try
                    %                         c_m=c_m+1;
                    % [r,lags]=xcorr(wells_corr(temp_low(m):temp_high(m)),tunnel_corr{k}(temp_low(m):temp_high(m)),'normalized');
                    [r,lags]=xcorr(wells_corr,tunnel_corr{k},'normalized');
%                     [temp_r(c_k), I]=max(r);
%                     temp_l(c_k)=lags(I);
                    %                          figure
                    %                         stem(lags,r);
                    %end
                    %catch
                    %    temp_l(m)=nan;
                    %    temp_r(m)=0;
                    %end
                %end
                [corr_r(c_k),I]=max(r);
                temp_lag(c_k)=lags(I);
                temp_cor_std(c_k)=std(temp_r);
                temp_lag_std(c_k)=std(temp_l);
                temp_n(c_k)=numel(tab.pl{1, fi}{1, k+5});

                names(c_k)={strcat(dir_lis(j)," & ",convertStringsToChars(matching_table{k+5,2}))};
                temp_fd(c_k)={convertCharsToStrings(Index_fd{fi,k+5})};

                switch Index_fd{fi,k+5}
                    case 'ff'
                        color=[color,"r"];
                    case 'fb'
                        color=[color,"b"];
                    case 'ff&fb'
                        color=[color,"g"];%"m"
                    case 'D'
                        color=[color,"g"];
                end
            end
            end
            if ~isempty(corr_r)
                c=c+1;
        [Bestcorr(c),I]=max(corr_r);
        Bestlagg(c)=temp_lag(I);
        all_name(c)={names};
        best_name(c)=names(I);
        fd_dr_best(c)={Index_fd{fi,I+5}};
        fd_de(c)={temp_fd};
        all_lag(c)={temp_lag};
        all_cor(c)={corr_r};
        all_l_std(c)={temp_lag_std};
        all_c_std(c)={temp_cor_std};
        all_n(c)={temp_n};
            end
    end
end
cor_tab.corr(fi)={Bestcorr};
cor_tab.lag(fi)={Bestlagg};
cor_tab.best_names(fi)={best_name};
cor_tab.all_names(fi)={all_name};
cor_tab.fd(fi)={fd_dr};
cor_tab.fd_best(fi)={fd_dr_best};
cor_tab.all_corr(fi)={all_cor};
cor_tab.all_corr_std(fi)={all_c_std};
cor_tab.all_lags_std(fi)={all_l_std};
cor_tab.all_lags(fi)={all_lag};
cor_tab.colors(fi)={color};
cor_tab.n(fi)={all_n};










end

%% Graph DG-CA3
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesTickLength',[0.04,0.08])
set(0,'defaultaxeslinewidth',2)
folders=["ECDGCA3CA1 19914 160127 160217 d21",...
    "ECDGCA3CA1 19914 160127 160303 d37",...
    "ECDGCA3CA1 19908 150729 150823 d25",...
    "ECDGCA3CA1 24088 160127 160302 d36",...
    "ECDGCA3CA1 19908 160518 160610 d22",...
    "ECDGCA3CA1 24574 160127 160303 d37",...
    "ECDGCA3CA1 19911 160518 160610 d22",...
    "ECDGCA3CA1 24574 160727 160818 d22",...
    "ECDGCA3CA1 19914 150805 150828 d25"];

cd 'graph'

if ~exist('CA3', 'dir')
    mkdir('CA3');
end

cd 'CA3'

for fi=1:9
    temp_names=[];
    temp_cor=[];
    temp_lag=[];
    temp_cor_std=[];
    temp_lag_std=[];
    temp_cn=[];

    for i=1:numel(cor_tab.all_names{1, fi})
        for j=1:numel(cor_tab.all_names{1, fi}{1,i})
            temp_names=[temp_names,convertCharsToStrings(cor_tab.all_names{1, fi}{1,i}{1,j})];
        end
        temp_cor=[temp_cor,cor_tab.all_corr{1, fi}{1,i}];
        temp_lag=[temp_lag,cor_tab.all_lags{1, fi}{1,i}];
        temp_cor_std=[temp_cor_std,cor_tab.all_corr_std{1, fi}{1,i}];
        temp_lag_std=[temp_lag_std,cor_tab.all_lags_std{1, fi}{1,i}];
        temp_cn=[temp_cn,cor_tab.n{1, fi}{1,i}];
    end
    figure('units','normalized','outerposition',[0 0 1 1])
    Tresh_cor=mean(temp_cor);
    c_s=0;
    c_p=0;
    for i=1:numel(temp_cor)
        if  cor_tab.colors{1, fi}(i)~="g" %& temp_cor(i)>=Tresh_cor
            h1=bar(categorical(temp_names(i)),temp_cor(i));
            set(h1, 'FaceColor', cor_tab.colors{1, fi}(i))
            hold on;
            errorbar(categorical(temp_names(i)),temp_cor(i),temp_cor_std(i)./sqrt(temp_cn(i)),'LineWidth',2);
            c_p=c_p+1;

        end
        if i<numel(temp_cor)
            temp_well_name1=convertStringsToChars(temp_names(i));
            temp_well_name2=convertStringsToChars(temp_names(i+1));
            if  length(temp_well_name1(1:strfind(temp_well_name1,'&'))) ~= length(temp_well_name2(1:strfind(temp_well_name2,'&'))) & c_p~=0
                c_s=c_s+1;
                c_p=0;
                h1=bar(categorical(string(c_s)),0);
                hold on;
            else
                if convertCharsToStrings(temp_well_name1(1:strfind(temp_well_name1,'&'))) ~= convertCharsToStrings(temp_well_name2(1:strfind(temp_well_name2,'&'))) & c_p~=0
                    c_s=c_s+1;
                    c_p=0;
                    h1=bar(categorical(string(c_s)),0);
                    hold on;
                end
            end
        elseif c_p~=0
            c_s=c_s+1;
            h1=bar(categorical(string(c_s)),0);
        end
    end
    ylabel('Correlations')
    ylim([0,1])
    %title(strcat("cor DG-CA3 ",folders(fi)))
    hold off;
    saveas(gcf,strcat("cor_DG_CA3 ",folders(fi)),'png')

    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:numel(temp_cor)
        if  cor_tab.colors{1, fi}(i)~="g" %& temp_cor(i)>=Tresh_cor % & temp_cor(i)>=Tresh_cor
            h1=bar(categorical(temp_names(i)),temp_lag(i));
            set(h1, 'FaceColor', cor_tab.colors{1, fi}(i))
            hold on;
        end

    end
    ylabel('lags(ms)')
%    title(strcat("lag DG-CA3 ",folders(fi)))
    hold off;
    saveas(gcf,strcat("lag_DG_CA3 ",folders(fi)),'png')
end


cd ..\..
%% data Cal for CA3-CA1
set(0,'defaultAxesFontSize',32)
dir_lis=["H8","H9","H10","H11","H12","J8","J9","J10","J11","J12","K8","K9","K10","K11","L8","L9","L10","M8","M9"];

load('test 32 CA3 400ms thr.mat')
wells=tab.rspin;


load('test 31 400ms thr.mat')
tunnel=tab.rspin;


is_cw=[1,1,0,1,0,1,0,1,0];
load('Index_fd_CA3_V2.mat');
for fi=1:9
    wells_temp=wells{1,fi};
    tunnel_temp=tunnel{1,fi};
    for i=6:15
        tunnel_corr(i-5)={tunnel_temp{1,i}};
    end
    if is_cw(fi)==1
        load matching_table_cw.mat
    else
        load matching_table_ccw.mat
    end
    c=0;
    Bestcorr=[];
    Bestlagg=[];
    fd_dr={};
    fd_dr_best={};
    all_lag={};
    all_cor={};
    all_l_std={};
    all_c_std={};
    all_n={};
    all_name={};
    best_name={};
    color=[];

    for j=1:19
        if ~isempty(wells_temp{1,j})
            
            wells_corr=wells_temp{1,j};
            corr_r=[];
            temp_lag=[];
            names={};
            temp_fd={};
            temp_lag_std=[];
            temp_cor_std=[];
            temp_n=[];
            c_k=0;

            for k=5:10 %1:5
                if ~isempty(tunnel_corr{k})
                    c_k=c_k+1;
                    temp_r=[];
                    temp_l=[];
                    c_m=0;
                    %                     for m=1:numel(tab.ph{1, fi}{1, k+5})
                    %                         temp_low=floor(tab.pl{1, fi}{1, k+5}./25);
                    %                         temp_high=floor(tab.ph{1, fi}{1, k+5}./25);
                    %                         temp_low(temp_low==0)=1;
                    %                         if temp_high(m)-temp_low(m)>=10
                    %                         %try
                    %                         c_m=c_m+1;
                    % [r,lags]=xcorr(wells_corr(temp_low(m):temp_high(m)),tunnel_corr{k}(temp_low(m):temp_high(m)),'normalized');
                    [r,lags]=xcorr(wells_corr,tunnel_corr{k},'normalized');
%                     [temp_r(c_k), I]=max(r);
%                     temp_l(c_k)=lags(I);
                    %                          figure
                    %                         stem(lags,r);
                    %end
                    %catch
                    %    temp_l(m)=nan;
                    %    temp_r(m)=0;
                    %end
                %end
                [corr_r(c_k),I]=max(r);
                temp_lag(c_k)=lags(I);
                temp_cor_std(c_k)=std(temp_r);
                temp_lag_std(c_k)=std(temp_l);
                temp_n(c_k)=numel(tab.pl{1, fi}{1, k+5});

                names(c_k)={strcat(dir_lis(j)," & ",convertStringsToChars(matching_table{k+5,2}))};
                temp_fd(c_k)={convertCharsToStrings(Index_fd{fi,k+5})};

                switch Index_fd{fi,k+5}
                    case 'ff'
                        color=[color,"r"];
                    case 'fb'
                        color=[color,"b"];
                    case 'ff&fb'
                        color=[color,"g"];%"m"
                    case 'D'
                        color=[color,"g"];
                end
            end
            end
            if ~isempty(corr_r)
                c=c+1;
        [Bestcorr(c),I]=max(corr_r);
        Bestlagg(c)=temp_lag(I);
        all_name(c)={names};
        best_name(c)=names(I);
        fd_dr_best(c)={Index_fd{fi,I+5}};
        fd_de(c)={temp_fd};
        all_lag(c)={temp_lag};
        all_cor(c)={corr_r};
        all_l_std(c)={temp_lag_std};
        all_c_std(c)={temp_cor_std};
        all_n(c)={temp_n};
            end
    end
end
cor_tab.corr(fi)={Bestcorr};
cor_tab.lag(fi)={Bestlagg};
cor_tab.best_names(fi)={best_name};
cor_tab.all_names(fi)={all_name};
cor_tab.fd(fi)={fd_dr};
cor_tab.fd_best(fi)={fd_dr_best};
cor_tab.all_corr(fi)={all_cor};
cor_tab.all_corr_std(fi)={all_c_std};
cor_tab.all_lags_std(fi)={all_l_std};
cor_tab.all_lags(fi)={all_lag};
cor_tab.colors(fi)={color};
cor_tab.n(fi)={all_n};










end
%% Graph CA3-CA1
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesTickLength',[0.02,0.04])
set(0,'defaultaxeslinewidth',2)
folders=["ECDGCA3CA1 19914 160127 160217 d21",...
    "ECDGCA3CA1 19914 160127 160303 d37",...
    "ECDGCA3CA1 19908 150729 150823 d25",...
    "ECDGCA3CA1 24088 160127 160302 d36",...
    "ECDGCA3CA1 19908 160518 160610 d22",...
    "ECDGCA3CA1 24574 160127 160303 d37",...
    "ECDGCA3CA1 19911 160518 160610 d22",...
    "ECDGCA3CA1 24574 160727 160818 d22",...
    "ECDGCA3CA1 19914 150805 150828 d25"];

cd 'graph'

if ~exist('CA3', 'dir')
    mkdir('CA3');
end

cd 'CA3'

for fi=1:9
    temp_names=[];
    temp_cor=[];
    temp_lag=[];
    temp_cor_std=[];
    temp_lag_std=[];
    temp_cn=[];

    for i=1:numel(cor_tab.all_names{1, fi})
        for j=1:numel(cor_tab.all_names{1, fi}{1,i})
            temp_names=[temp_names,convertCharsToStrings(cor_tab.all_names{1, fi}{1,i}{1,j})];
        end
        temp_cor=[temp_cor,cor_tab.all_corr{1, fi}{1,i}];
        temp_lag=[temp_lag,cor_tab.all_lags{1, fi}{1,i}];
        temp_cor_std=[temp_cor_std,cor_tab.all_corr_std{1, fi}{1,i}];
        temp_lag_std=[temp_lag_std,cor_tab.all_lags_std{1, fi}{1,i}];
        temp_cn=[temp_cn,cor_tab.n{1, fi}{1,i}];
    end
    figure('units','normalized','outerposition',[0 0 1 1])
    Tresh_cor=mean(temp_cor);
    c_s=0;
    c_p=0;
    for i=1:numel(temp_cor)
        if  cor_tab.colors{1, fi}(i)~="g" %& temp_cor(i)>=Tresh_cor
            h1=bar(categorical(temp_names(i)),temp_cor(i));
            set(h1, 'FaceColor', cor_tab.colors{1, fi}(i))
            hold on;
            errorbar(categorical(temp_names(i)),temp_cor(i),temp_cor_std(i)./sqrt(temp_cn(i)),'LineWidth',2);
            c_p=c_p+1;
            er.Color = [0 0 0];
            er.LineStyle = 'none';



        end
        if i<numel(temp_cor)
            temp_well_name1=convertStringsToChars(temp_names(i));
            temp_well_name2=convertStringsToChars(temp_names(i+1));
            if  length(temp_well_name1(1:strfind(temp_well_name1,'&'))) ~= length(temp_well_name2(1:strfind(temp_well_name2,'&'))) & c_p~=0
                c_s=c_s+1;
                c_p=0;
                h1=bar(categorical(string(c_s)),0);
                hold on;
            else
                if convertCharsToStrings(temp_well_name1(1:strfind(temp_well_name1,'&'))) ~= convertCharsToStrings(temp_well_name2(1:strfind(temp_well_name2,'&'))) & c_p~=0
                    c_s=c_s+1;
                    c_p=0;
                    h1=bar(categorical(string(c_s)),0);
                    hold on;
                end
            end
        elseif c_p~=0
            c_s=c_s+1;
            h1=bar(categorical(string(c_s)),0);
        end

    end
    ylabel('Correlations')
    ylim([0,1])
    title(strcat("cor CA3-CA1 ",folders(fi)))
    hold off;
    saveas(gcf,strcat("cor_CA3_CA1 ",folders(fi)),'png')

    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:numel(temp_cor)
        if  cor_tab.colors{1, fi}(i)~="g" %& temp_cor(i)>=Tresh_cor %
            h1=bar(categorical(temp_names(i)),temp_lag(i));
            set(h1, 'FaceColor', cor_tab.colors{1, fi}(i))
            hold on;
            errorbar(categorical(temp_names(i)),temp_lag(i),temp_lag_std(i)./sqrt(temp_cn(i)),'LineWidth',2);
            c_p=c_p+1;
            er.Color = [0 0 0];
            er.LineStyle = 'none';

        end

        if i<numel(temp_cor)
            temp_well_name1=convertStringsToChars(temp_names(i));
            temp_well_name2=convertStringsToChars(temp_names(i+1));
            if  length(temp_well_name1(1:strfind(temp_well_name1,'&'))) ~= length(temp_well_name2(1:strfind(temp_well_name2,'&'))) & c_p~=0
                c_s=c_s+1;
                c_p=0;
                h1=bar(categorical(string(c_s)),0);
                hold on;
            else
                if convertCharsToStrings(temp_well_name1(1:strfind(temp_well_name1,'&'))) ~= convertCharsToStrings(temp_well_name2(1:strfind(temp_well_name2,'&'))) & c_p~=0
                    c_s=c_s+1;
                    c_p=0;
                    h1=bar(categorical(string(c_s)),0);
                    hold on;
                end
            end
        elseif c_p~=0
            c_s=c_s+1;
            h1=bar(categorical(string(c_s)),0);
        end


    end
    ylabel('lags(ms)')
    title(strcat("lag CA3-CA1 ",folders(fi)))
    hold off;
    saveas(gcf,strcat("lag_CA3_CA1 ",folders(fi)),'png')
end


cd ..\..





