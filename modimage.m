function bdctimgs=modimage(spimgs,targetimgs)
%function that modifies image

bdctimgs=zeros(size(spimgs));
T=10;
c=clock;
filename=strcat('tempresult',int2str(c(2)),int2str(c(3)),int2str(c(4)),int2str(c(5)));
logfilename=strcat('logfile',int2str(c(2)),int2str(c(3)),int2str(c(4)),int2str(c(5)));
logfileid=fopen(logfilename,'a');
for i=1:size(spimgs,1)
    spimg=spimgs(i,:);
    spimg=reshape(spimg,128,128);
    targetimg=targetimgs(i,:);
    targetimg=reshape(targetimg,128,128);
    bdctimg=hopfieldnet(spimg,targetimg,T);
    bdctimgs(i,:)=bdctimg(:)';
    save(filename,bdctimgs);
    fprintf(logfileid,'%s: finished processing image %d\n',strcat(int2str(c(2)),int2str(c(3)),int2str(c(4)),int2str(c(5))),i);
end

fclose(logfileid);

