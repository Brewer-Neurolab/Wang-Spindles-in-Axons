%% FF
Fs = 1000;            % Sampling frequency 1000
T = 1/Fs;             % Sampling period
L = 300000;            % Length of signal 300000
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;
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
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];
    temp_setff=[];temp_setfb=[];
    c_ff=0;
    for fi=1:9
        
        tmp=tab.rspin{1, fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff' & ~isempty(tmp{k+(regi-1)*5})
                Y = fft(tmp{k+(regi-1)*5});
                P2 = abs(Y/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                P1=P1(f<=16 & f>=10);
                c_ff=c_ff+1;
                temp_setff(c_ff,:)=P1;

            elseif Index_fd{fi,k+(regi-1)*5}=='fb' & ~isempty(tmp{k+(regi-1)*5})
                Y = fft(resample(tmp{k+(regi-1)*5},1,1));
                P2 = abs(Y/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                
                P1=P1(f<=16 & f>=10);

            end
        end
 
    end

    subplot(2,2,figOrder(regi))
    if ~isempty(temp_setff)
    for i=1:numel(P1)
    setff(i)=mean(temp_setff(:,i));
    end
    end
    if regi~=2
    plot(10:6/(numel(setff)/100):16,resample(setff,1,100), 'LineWidth',2)
    else
    plot([],[], 'LineWidth',2)
    end
    xticks([10:16])
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    yticks([0:0.005:0.025])
    ylim([0,0.025])
    title(well_regList(regi))
    set(gca,'fontsize',18)





end

saveas(gcf,strcat("FFT_062523",'ff'),'png')

cd ..
%% FB
T = 1/Fs;             % Sampling period
L = 300000;             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;
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
co=0;
for regi = 1:4
    coliff=1; histcountff=[];
    colifb=1; histcountfb=[];
    temp_setff=[];temp_setfb=[];
    c_fb=0;
    for fi=1:9
        
        tmp=tab.rspin{1, fi};
        for k=1:5
            if Index_fd{fi,k+(regi-1)*5}=='ff' & ~isempty(tmp{k+(regi-1)*5})
                Y = fft(tmp{k+(regi-1)*5});
                P2 = abs(Y/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                P1=P1(f<=16 & f>=10);

            elseif Index_fd{fi,k+(regi-1)*5}=='fb' & ~isempty(tmp{k+(regi-1)*5})
                                Y = fft(tmp{k+(regi-1)*5});
                P2 = abs(Y/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                
                P1=P1(f<=16 & f>=10);
                c_fb=c_fb+1;
                temp_setfb(c_fb,:)=P1;


            end
        end
 
    end

    for i=1:numel(P1)
    setfb(i)=mean(temp_setfb(:,i));
    end

    subplot(2,2,figOrder(regi))
    plot(10:6/(numel(setfb)/100):16,resample(setfb,1,100),'LineWidth',2)
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
  xticks([10:16])
  ylim([0,0.025])
     yticks([0:0.005:0.025])
    title(well_regList(regi))
    set(gca,'fontsize',18)





end

saveas(gcf,strcat("FFT_062523",'fb'),'png')

cd ..