 clear all; %close all; clc
load matching_table_ccw.mat
folders=["ECDGCA3CA1 19914 160127 160217 d21 5minspont000_mat_files",...
    "ECDGCA3CA1 19914 160127 160303 d37 5minspont000_mat_files",...
"ECDGCA3CA1 19908 150729 150823 d25 5minspont000_mat_files",...
"ECDGCA3CA1 24088 160127 160302 d36 5minspont000_mat_files",...
"ECDGCA3CA1 19908 160518 160610 d22 5minspont000_mat_files",...  
"ECDGCA3CA1 24574 160127 160303 d37 5minspont000_mat_files",...  
"ECDGCA3CA1 19911 160518 160610 d22 5minspont000_mat_files",...
"ECDGCA3CA1 24574 160727 160818 d22 5minspont000_mat_files",...  
"ECDGCA3CA1 19914 150805 150828 d25 5minspont000_mat_files"];
s_folders=["ECDGCA3CA1 19914 160127 160217 d21 5minspont0001_mat_files",...
    "ECDGCA3CA1 19914 160127 160303 d37 5minspont0001_mat_files",...
"ECDGCA3CA1 19908 150729 150823 d25 5minspont0001_mat_files",...
"ECDGCA3CA1 24088 160127 160302 d36 5minspont0001_mat_files",...
"ECDGCA3CA1 19908 160518 160610 d22 5minspont0001_mat_files",...  
"ECDGCA3CA1 24574 160127 160303 d37 5minspont0001_mat_files",...  
"ECDGCA3CA1 19911 160518 160610 d22 5minspont0001_mat_files",...
"ECDGCA3CA1 24574 160727 160818 d22 5minspont0001_mat_files",...  
"ECDGCA3CA1 19914 150805 150828 d25 5minspont0001_mat_files"];
Chanel_Name=[];
tab.Fol={};
tab.In={};
tab.Len={};
tab.IC={};
is_cw=[1,1,0,1,0,1,0,1,0];
qcw=[];
tic
for j=5:9 %(1:9)
    if is_cw(j)==1
        load matching_table_cw.mat
        qcw=[qcw,'T'];
    else
        load matching_table_ccw.mat
        qcw=[qcw,'F'];
    end
    Interval(1:20)=cell(1,20);
    Length(1:20)=cell(1,20);
    Integral_conv(1:20)=cell(1,20);
    for i=6:20 %1:20
        %load Spindle-
        cd('Spindle')
        cd(folders(j))
        temp_dir= convertStringsToChars(matching_table{i,2})
        temp_dir(strfind(temp_dir,'-'):end)=[];
        %temp_dir(1:strfind(temp_dir,'-'))=[];
        Chanel_Name{i}={temp_dir};
        temp_dir1=[temp_dir,'.mat'];
        data=load (temp_dir1);
        data1=data.data;
        cd ..\..
        %load Spike
        cd('Spike')
        cd(s_folders(j))
        temp_dir= convertStringsToChars(matching_table{i,2});
        temp_dir(strfind(temp_dir,'-'):end)=[];
        Chanel_Name{i}={temp_dir};
        temp_dir2=strcat('times_',temp_dir,'.mat');
        temp_dir=[temp_dir,'.mat'];
        data=load (temp_dir);
        data2=data.data;

        try
 data=load (temp_dir2, 'cluster_class');
        data3=data.cluster_class  (:,2);
        data3=data3(data.cluster_class(:,1)==1);
        catch 
            data3=[];
        end
        cd ..\..

        %calculation

        [In,Len,IC,pl,ph,pc,re_spike,re_spindle,peak,s_angle,spike_count]=Spindle(data1,data2,data3);
        Interval(i)={In};
        Length(i)={Len};
        Integral_conv(i)={IC};
        po_l(i)={pl};
        po_h(i)={ph};
        pcs(i)={pc};
        datak(i)={data3};
        re_spik(i)={re_spike};
        re_spin(i)={re_spindle};
        peak_amp(i)={peak};
        an(i)={s_angle};
        sc(i)={spike_count};
        %         conv(i)={conv_s};
        %         sig(i)={s_out};
    end
    tab.Fol(j)={folders(j)};
    tab.In(j)={Interval};
    tab.Len(j)={Length};
    tab.IC(j)={Integral_conv};
    tab.pl(j)={po_l};
    tab.ph(j)={po_h};
    tab.pc(j)={pcs};
    tab.data3(j)={datak};
    tab.rspike(j)={re_spik};
    tab.rspin(j)={re_spin};
    tab.peak_amp(j)={peak_amp};
    tab.angle(j)={an};
    tab.spik_n(j)={sc};
end
