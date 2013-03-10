function bdctimgs=modimage2(spimgs,targetimgs)
%function that modifies image, parfor version

bdctimgs=zeros(size(spimgs));
T=10;
c=clock;
filename=strcat('tempresult',int2str(c(2)),int2str(c(3)),int2str(c(4)),int2str(c(5)));
logfilename=strcat('logfile',int2str(c(2)),int2str(c(3)),int2str(c(4)),int2str(c(5)));
M=sqrt(size(spimgs,2));
for i=1:size(spimgs,1)
    spimg=spimgs(i,:);
    spimg=reshape(spimg,M,M);
    targetimg=targetimgs(i,:);
    targetimg=reshape(targetimg,M,M);
    bdctimg=hopfieldnet(spimg,targetimg,T);
    bdctimgs(i,:)=bdctimg(:)';
    save(filename,'bdctimgs');
    logfileid=fopen(logfilename,'a');
    c=clock;
    fprintf(logfileid,'%s: finished processing image %d\n',strcat(int2str(c(2)),int2str(c(3)),int2str(c(4)),int2str(c(5))),i);
    fclose(logfileid);
end



