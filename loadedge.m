function edgeimg=loadedge(imgpath,imgname)
%load edge image given an image

edgeimg=imread([imgpath filesep 'edgemask' filesep imgname(1:end-4) '_edgemask.jpg']);

