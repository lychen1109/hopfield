function [aupicset1,aupicset2,sppicset1,sppicset2]=filedivide
%filedivide divide files into 4 category

pics={'canong3_02','canong3_05','canong3_08',...
    'nikond70_02','nikond70_05','nikond70_08','nikond70_11',...
    'canonxt_02','canonxt_05','canonxt_08','canonxt_11','canonxt_14','canonxt_17','canonxt_20','canonxt_23','canonxt_26','canonxt_29','canonxt_32','canonxt_35','canonxt_38'};
kodakpics={'kodakdcs330_01_sub_01.tif','kodakdcs330_02_sub_01.tif','kodakdcs330_03_sub_01.tif'};

picidx=randperm(20);
aupicset1=cell(90,1);
aupicset2=cell(93,1);

%blocks divide for the first 10 pictures
for i=1:10
    blockidx=randperm(9);
    aupicset1{(i-1)*4+1}=[pics{picidx(i)} '_sub_0' int2str(blockidx(1)) '.tif'];
    aupicset1{(i-1)*4+2}=[pics{picidx(i)} '_sub_0' int2str(blockidx(2)) '.tif'];
    aupicset1{(i-1)*4+3}=[pics{picidx(i)} '_sub_0' int2str(blockidx(3)) '.tif'];
    aupicset1{(i-1)*4+4}=[pics{picidx(i)} '_sub_0' int2str(blockidx(4)) '.tif'];
    aupicset2{(i-1)*5+1}=[pics{picidx(i)} '_sub_0' int2str(blockidx(5)) '.tif'];
    aupicset2{(i-1)*5+2}=[pics{picidx(i)} '_sub_0' int2str(blockidx(6)) '.tif'];
    aupicset2{(i-1)*5+3}=[pics{picidx(i)} '_sub_0' int2str(blockidx(7)) '.tif'];
    aupicset2{(i-1)*5+4}=[pics{picidx(i)} '_sub_0' int2str(blockidx(8)) '.tif'];
    aupicset2{(i-1)*5+5}=[pics{picidx(i)} '_sub_0' int2str(blockidx(9)) '.tif'];
end

%blocks divide for the next 10 pictures
for i=11:20
    blockidx=randperm(9);
    aupicset1{40+(i-11)*5+1}=[pics{picidx(i)} '_sub_0' int2str(blockidx(1)) '.tif'];
    aupicset1{40+(i-11)*5+2}=[pics{picidx(i)} '_sub_0' int2str(blockidx(2)) '.tif'];
    aupicset1{40+(i-11)*5+3}=[pics{picidx(i)} '_sub_0' int2str(blockidx(3)) '.tif'];
    aupicset1{40+(i-11)*5+4}=[pics{picidx(i)} '_sub_0' int2str(blockidx(4)) '.tif'];
    aupicset1{40+(i-11)*5+5}=[pics{picidx(i)} '_sub_0' int2str(blockidx(5)) '.tif'];
    aupicset2{50+(i-11)*4+1}=[pics{picidx(i)} '_sub_0' int2str(blockidx(6)) '.tif'];
    aupicset2{50+(i-11)*4+2}=[pics{picidx(i)} '_sub_0' int2str(blockidx(7)) '.tif'];
    aupicset2{50+(i-11)*4+3}=[pics{picidx(i)} '_sub_0' int2str(blockidx(8)) '.tif'];
    aupicset2{50+(i-11)*4+4}=[pics{picidx(i)} '_sub_0' int2str(blockidx(9)) '.tif'];
end

%add last three kodak pictures into aupicset2

aupicset2(91:93)=kodakpics;

sppairs={'canong3_canonxt','canong3_kodakdcs330','canong3_nikond70','canonxt_kodakdcs330','nikond70_canonxt','nikond70_kodakdcs330'};
sppicset1=cell(90,1);
sppicset2=cell(90,1);
for i=1:length(sppairs)
    blockidx=randperm(30);
    for j=1:15
        if blockidx(j)<10
            sppicset1{(i-1)*15+j}=[sppairs{i} '_sub_0' int2str(blockidx(j)) '.tif'];            
        else
            sppicset1{(i-1)*15+j}=[sppairs{i} '_sub_' int2str(blockidx(j)) '.tif'];            
        end
        if blockidx(j+15)<10
            sppicset2{(i-1)*15+j}=[sppairs{i} '_sub_0' int2str(blockidx(j+15)) '.tif'];
        else
            sppicset2{(i-1)*15+j}=[sppairs{i} '_sub_' int2str(blockidx(j+15)) '.tif'];
        end
    end
end
        


    