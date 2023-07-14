clear;
close all;

% Visible image enhancement path
save_dir='.\envis\';
if exist(save_dir)==0 
    mkdir(save_dir);
else
    disp('dir is exist');
end

% Fusion path
save_dir1='.\fusion\';
if exist(save_dir1)==0 
    mkdir(save_dir1);
else
    disp('dir is exist');
end

for i = 2:2 % TNO dataset
    PathIr          = [ 'IV_images\IR' ,        num2str(i) ,       '.png' ];
    PathVIS          = [ 'IV_images\VIS' ,        num2str(i) ,       '.png' ];
    %% 1. Visible image enhancement
    C = imread(PathVIS);
    C = im2double(C);
    envis = Enhancement(C);
    imwrite( envis ,[save_dir ['envis_',num2str(i),'.png']]);
%% 2.1 Prepare the infrared image and the enhanced visible image
    A  = imread(PathIr);
    A = im2double(A);   
    A = wmcFilter(A,3);
    B = im2double(envis);
%     B = wmcFilter(B,3);
%% 2.2 Put infrared and visible images into a sequence
    imgSeqColor = zeros(size(A,1),size(A,2),3,2); 
    if size(A,3)==1
    imgSeqColor(:,:,:,1) = cat(3,A,A,A);
    else
    imgSeqColor(:,:,:,1) = A;
    end

    if size(B,3)==1
    imgSeqColor(:,:,:,2) = cat(3,B,B,B);
    else
    imgSeqColor(:,:,:,2) = B;
    end
    tic
    [r,c,ch,n] = size(imgSeqColor); 
   
    %% 3.Multiscale processing
    %% the fine  scale
    r1=4;
    [ Detail,i_mean1,aa1,N1] = scale_fine_denoise(imgSeqColor,r1); 
    %% the intermediate  scale
    [w,h,~,~]=size(imgSeqColor);
    
    nlev = floor(log(min(w,h)) / log(2))-5; 
   
    D2 = cell(nlev,1);
    aa2= cell(nlev,1);
    N2= cell(nlev,1);
   
    r2=4;
        for ii=1:nlev
            [ D2{ii},i_mean2,aa2{ii},N2{ii}] = scale_interm(i_mean1,r2);
            i_mean1=i_mean2;
        end


    %% the coarsest  scale
    r3=4;

    [fI3,i_mean3,aa3,N3] = zyq_scale_coarse(i_mean2,r3);

    %% reconstruct
    %% Intermediate layers
        for ii=nlev:-1:1
            temp=aa2{ii};
            fI=zeros(size(temp));
            fI(1:2:size(temp,1),1:2:size(temp,2))=fI3;
            B2=boxfilter(fI, r2)./ N2{ii}+D2{ii};
            fI3=B2;
        end
    %% finest layers
    fI=zeros(size(aa1));
    fI(1:2:size(aa1,1),1:2:size(aa1,2))=B2;
    B1=boxfilter(fI, r1)./ N1;  
    C_out=repmat(B1,[1 1 3])+Detail; 
    toc
    figure,imshow(C_out),title('fusion');
%     imwrite( C_out ,[save_dir1 ['F_',num2str(i),'.png']]);
%     end
end