function histhacksp_test(pathau,pathsp,aufilename,spfilename,T)
%test funciton of histhacksp

auimg=imread([pathau filesep aufilename]);
ycbcr=rgb2ycbcr(auimg);
cbau=double(ycbcr(:,:,2));

spimg=imread([pathsp filesep spfilename]);
ycbcr=rgb2ycbcr(spimg);
cbsp=double(ycbcr(:,:,2));

tpmopt=1;
tmtarget=tpm1(cbau,T,tpmopt);
tmori=tpm1(cbsp,T,tpmopt);
fprintf('Original distance between features %g\n',norm(tmori(:)-tmtarget(:)));

cbnew=histhacksp(cbsp,cbau,4,T,[]);
tm=tpm1(cbnew,T,tpmopt);
fprintf('Final distance between features %g\n',norm(tm(:)-tmtarget(:)));