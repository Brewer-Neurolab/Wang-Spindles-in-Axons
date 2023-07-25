function [gap,len,integral_conv,pl,ph,pc,re_spike,re_spindle,peak,s_angle,spike_count]=Spindle(data,data2,data3)
% out put zeroing for noise channel
fc1=10;%10
fc2=16;%16
fs=25000;
spindle_only=[];
gap=[];
len=[];
spike_count=[];
conv_s=[];
integral_conv=[];
s_out=[];
re_spike=[];
re_spindle=[];
peak=[];
s_angle=[];
pc=[];
pl=[];
ph=[];
T_spindle=[];
T_spike=[];
spike_array(1:7500000)=0;
s_in=data;
% [A,B,C,D] = butter(6,5/(fs/2),'high');
% %[B,A] = butter(3,[fc1/(fs/2),fc2/(fs/2)],'bandpass');
% [sos,g] = ss2sos(A,B,C,D);
% s_out1=filtfilt(sos,g,s_in);
% [A,B,C,D] = butter(8,100/(fs/2),'low');
% [sos,g] = ss2sos(A,B,C,D);
% s_out=filtfilt(sos,g,s_out1);
% s_in=medfilt1(s_out,10);
% signal frequency fitlerting
 % sampling frequency
x=[1/fs:1/fs:300]; % time axis
N=0.5*fs;
% [A,B,C,D] = butter(3,100/(fs/2),'low');
% [B,A] = butter(3,[fc1/(fs/2),fc2/(fs/2)],'bandpass');
% [sos,g] = ss2sos(A,B,C,D);
%s_out2=filtfilt(sos,g,s_in);
%% Fitler for spindle
[A,B,C,D] = butter(3,[1/(fs/2),50/(fs/2)],'bandpass');
[sos,g] = ss2sos(A,B,C,D);
s_625=filtfilt(sos,g,s_in);



%[A,B,C,D] = butter(3,[fc1/(fs/2),fc2/(fs/2)],'bandpass'); % 3
[A,B,C,D] = butter(6,fc1/(fs/2),'high');
%[B,A] = butter(3,[fc1/(fs/2),fc2/(fs/2)],'bandpass');
[sos,g] = ss2sos(A,B,C,D);
s_out1=filtfilt(sos,g,s_in);
[A,B,C,D] = butter(8,fc2/(fs/2),'low');
[sos,g] = ss2sos(A,B,C,D);
s_out=filtfilt(sos,g,s_out1);

fsamp = fs;
fcuts = [0.1 10 16 25000];
mags = [0 1 0];
devs = [0.01 0.05 0.01];

% [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
% n = n + rem(n,2);
% hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
% s_out=filtfilt(hh,1,s_in);

% [A,B,C,D] = butter(6,[fc1/(fs/2),fc2/(fs/2)],'bandpass'); % 3
% [sos,g] = ss2sos(A,B,C,D);
% s_out=filtfilt(sos,g,s_out2);
% s_out=sosfilt(sos,s_in);
% sf=fft(s_out);
s_out=detrend(s_out);
% envelope generation
h_s=hilbert(s_out);
h_s=abs(h_s);
gw=gausswin(numel(x)/2500,2.5); % questionable convolution
%conv_s=conv(h_s,gw,'same')/(numel(x)/6250);
conv_s=h_s;
% spindle detection
S_d=std(s_out);
me=mean(s_out);

upper_thr=me+2.5*S_d;
lower_thr=me+1.5*S_d;
% Once excedes the Upper threshold, classify as spindle, and after excedes the lower threshold for at least 500ms, classify as spindles
abs_lowethr=2; %uV %10,2
I_l=find(conv_s>lower_thr & conv_s>abs_lowethr);
I_h=find(conv_s>upper_thr & conv_s>abs_lowethr);
c=1;
j=1;
I_l=[0,I_l];
p_s=[];
p_sh=[];
% detected region length detection for lower thresh hold
for i=2:(numel(I_l)-2)
    if I_l(i+1)==I_l(i)+1
        if I_l(i)~=I_l(i-1)+1
            if j>=12500
                c=c+1;
            elseif ~isempty(p_s)
                p_s(c,:)=[];
            end
            j=1;
        end
        p_s(c,j)=I_l(i);
        j=j+1;
    end
end





spindel=[];
t=1;
he_l=[];
po_l=[];
% see if noise channle
if isempty(p_s)
    figure
    plot(x,s_out)
    xlabel('time(s)')
    ylabel('Amplitude(uV)')
    hold on;
    plot(x,conv_s);
%    pause
    return
else
    % spindle recording
    for i=1:c
        if sum(p_s(i,:)~=0)>=N
            temp=p_s(i,:);
            temp=temp(temp~=0);
            he_l(t)=min(conv_s(temp));
            po_l(t)=min(temp);
            t=t+1;
            spindel=[spindel,temp];
        end
    end
end

t=t-1;
d=1;
po_h(1)=I_h(1);
sum_h=[];
num(1)=1;
sum_h(1)=conv_s(I_h(1));
min_h(1)=conv_s(I_h(1));
% over lap for high treshold spindle
for i=2:numel(I_h)
    if I_h(i)-I_h(1)==1
        sum_h(d)=sum_h(d)+conv_s(I_h(i));
        min_h(d)=min(min_h(d),conv_s(I_h(i)));
        num(d)=num(d)+1;
    else
        d=d+1;
        po_h(d)=I_h(i);
        sum_h(d)=conv_s(I_h(i));
        min_h(d)=conv_s(I_h(i));
        num(d)=1;
    end
end
he_h=min_h;
i=1;
% position detection
if isempty(po_l)
    pl_l=[];ph_l=[];
else
    [pl_l,ph_l]=fun_po(conv_s,he_l,po_l,t);
end
if isempty(po_h)
    pl_h=[];ph_h=[];
else
    [pl_h,ph_h]=fun_po(conv_s,he_h,po_h,d);
end
if isempty(I_h)
else
    spindel=[spindel,I_h];
end
spindel=spindel(spindel~=0);

%Inter spindle interval detection (Overlap)
pl=sort([pl_l,pl_h]);
ph=sort([ph_l,ph_h]);
re=1;
% disp(numel(pl))
mark=0;
while mark~=1
    gap=pl(2:end)-ph(1:end-1);
    gap=[8000,gap,8000];
    pl=pl(gap(1:end-1)>7500);
    ph=ph(gap(2:end)>7500);
    if min(gap>0)
        mark=1;
    end
end
pl=unique(pl);
ph=unique(ph);

len=ph-pl;
ph(len<=10000)=[];
pl(len<=10000)=[];
len(len<=10000)=[];
spindel_index=[1:1:numel(x)];
gap=pl(2:end)-ph(1:end-1);

if numel(pl)<=2
gap=[];len=[];integral_conv=[];pl=[];ph=[];pc=[];re_spike=[];
re_spindle=[];peak=[];s_angle=[];spike_count=[];
    return
end


%% Spike
index=floor(data3/1000*fs);
spike_array(index)=1;
%
% f_spike = bandpass(data2,[300,10000],fs);
% h_sk=abs(hilbert(f_spike));
% gw=gausswin(numel(x)/2500,2.5)/190;
 conv_spike=conv(spike_array,gw,'same');

a=s_out(index);
pc=numel(a(a>0))/numel(index);

%% find max in each spindle
index_ws=[];
for i=1:numel(pl)
    [peak(i),peak_po(i)]=max(conv_s(spindel_index>=pl(i) & spindel_index<=ph(i)));
    peak_po(i)=peak_po(i)+pl(i)-1;
    tmpS=index(index<=ph(i) & index>=pl(i));
    index_ws=[index_ws,tmpS'];
    spike_count(i)=numel(tmpS);
end




%% angle
te=hilbert(s_out);
sig=angle(te);
s_angle=sig(index_ws);





%% cal intergal
    spindle_only=conv_s;
    %spindle_only=s_out;
    spindle_only(1:pl(1))=0;
    pl(numel(pl)+1)=numel(conv_s);
for i=1:numel(ph)
    integral_conv(i)=sum(conv_s(pl(i):ph(i)));%trapz(x(pl(i):ph(i)),conv_1(pl(i):ph(i)));%
    %     temp=abs(s_out(pl(i):ph(i)));
    integral_raw(i)=sum(s_out(pl(i):ph(i)));%trapz(x(pl(i):ph(i)),temp);%
    %     T_spindle=[T_spindle,h_s(pl(i):ph(i))];
    %     T_spike=[T_spike,h_sk(pl(i):ph(i))];
    spindle_only(ph(i):pl(i+1))=0;
end

spindle_only=resample(spindle_only,1,25);
re_spindle=spindle_only;
% re_spindle=resample(s_out,1,25);
re_spike=resample(conv_spike,1,25);
%% ploting
%{
tiledlayout(3,1)

ax1 = nexttile;
plot(x,data)

ax2 = nexttile;
plot(x,s_in)

%ax2 = nexttile;
%plot(x,s_out2)
ax3 = nexttile;
plot(x,s_out)
hold on
plot(x,h_s)
xlabel('time(s)')
ylabel('Amplitude(uV)')
plot(x,conv_s);
% plot(spindel./fs,conv_s(spindel),'o')
plot(pl./fs,conv_s(pl),'o')
plot(ph./fs,conv_s(ph),'o')
plot(x(index_ws),s_out(index_ws),'o')
legend('10-16HZ','Hilibert Transform','envelope','start','end')
lgd.FontSize = 60;
%
linkaxes([ax1 ax2 ax3],'x')
%}




figure
set(0,'defaultAxesFontSize',32)
set(0,'defaultAxesTickLength',[0.01,0.02])
set(0,'defaultaxeslinewidth',2)

tiledlayout(4,1)
%plot 10-16hz
ax1 = nexttile;
plot(x,s_in,'LineWidth',2)
xlabel('time(s)')
ylabel('Amplitude(uV)')

%ax2 = nexttile;
%plot(x,s_out2)
ax2 = nexttile;
plot(x,s_out,'LineWidth',2)
hold on
plot(x,h_s,'LineWidth',2)
xlabel('time(s)')
ylabel('Amplitude(uV)')
% plot(x,conv_s);
% plot(spindel./fs,conv_s(spindel),'o')
plot(pl./fs,conv_s(pl),'o')
plot(ph./fs,conv_s(ph),'o')
ylim([-700,700])
%plot(x(index_ws),s_out(index_ws),'o')
%legend('10-16HZ','Hilibert Transformed envelope','start','end')
lgd.FontSize = 60;
ax3 = nexttile(3,[2,1]);
re_sig=resample(s_in,1,25);
[wt,f]=cwt(re_sig,1000,FrequencyLimits=[0 350]);
% 
% pcolor(1:numel(re_sig)/1000,log2(f1),abs(wt1));
% shading interp;
fig = plot_cwt(wt, f, 0:0.001:300, [], 0);
ylim([-2,8])
linkaxes([ax1 ax2 ax3],'x')

% yline(lower_thr);
% yline(upper_thr);
% stairs([0.001:0.001:300],re_spindle)
% Spikes
% ax3 = nexttile;
% 
% plot(x,spike_array)
% hold on;
% nu=[1:numel(x)];
% % plot(data3/1000,f_spike(index),'ro')
% plot(x,conv_spike)
% legend('Spike train','Envelope')
% 
% ax4 = nexttile;
% plot(x,rad2deg(sig))
% xlabel('time(s)')
% ylabel('angle(radius)')
% hold on;
% plot(x(index_ws),rad2deg(s_angle),'o')
% legend('Phase Angle','Spike position')
% 
% linkaxes([ax1 ax2 ax3 ax4],'x')

% stairs([0.001:0.001:300],re_spike)
%}
%}
end



%%Functions
function [ hist_object,hist_plot, hist_count, xCenter ] = log_binned__histogram( vec, binning_prarmeter, norm_method, if_plot )
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

function fig = plot_cwt(wt, wfreq, ts, cscale, norm)
% fig = plot_cwt(wt, wfreq, ts, if_log)
% Plot the output of CWT output 
% 
% Inputs -
% wt - wavelet coefficients 
% wfreq - freq vector
% ts - time series
% if_log - if magnitude in dB 
% 
% Created on Wed Oct 13 13:12:56 2021
% @author: yashvakilna

%% default parameters
if ~exist('if_log','var') 
    norm=0;
elseif isempty(norm)
    norm=0;
end

if ~exist('yscale','var')
    cscale = [1 100];
elseif isempty(cscale)
    cscale = [1 100];
end

if norm
    imagesc(ts, log2(wfreq), 10.*log10(abs(wt)))
else
    imagesc(ts, log2(wfreq), (abs(wt)))
end
%c=jet;
ax=gca;
ax.YDir = 'normal';
yticks = 2.^(round(log2(min(wfreq))):round(log2(max(wfreq))));
% yticks = [0.5,8,128];
hcol = colorbar; colormap(jet);
hcol.Location='east';

if norm
    hcol.Color='black';
    hcol.Label.String = "Magnitude";
    % set c-axis
    cah=caxis;
    caxis([cah(1)*0.05 cah(2)*0.95]);
else
%     hcol.Label.String = "dB/Hz";
    hcol.Label.String = "Wavelet power (uV^2)";
    hcol.Color='white';
    hcol.Label.FontSize=20;
end
set(ax,'YLim',log2(cscale),...
    'YTick',log2(yticks(:)), ...
    'YTickLabel',num2str(sprintf('%g\n',yticks)), ...
    'layer','top');

 ylabel('Frequency(Hz)')
     xlabel('Time (s)');
fig=gcf;

end
