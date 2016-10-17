function [lobarvalues, leftdicecoeff, rightdicecoeff] = lobarquantification( imagepathCT, imagepathQ, imagepathV, template, outputlabels)
% minor edits in comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set some global vbles
global displaydebugimages;
displaydebugimages=3;
startsliceCT=1;
useSDFimages=1;
lungmaxthreshold=100;
numerosions=0;

conn(:,:,1)=[0 0 0;0 1 0;0 0 0];
conn(:,:,2)=[0 1 0;1 1 1;0 1 0];
conn(:,:,3)=[0 0 0;0 1 0;0 0 0];
mincomponentsize=10000;
lefttranslateflag=0;
leftrigidflag=0;
leftsimilarityflag=0;
leftaffineflag=0;
leftnonlinearflag=1;
righttranslateflag=0;
rightrigidflag=0;
rightsimilarityflag=0;
rightaffineflag=0;
rightnonlinearflag=0;
lobarvalues=0;
leftdicecoeff=-1;
rightdicecoeff=-1;


% set output values in case early exit
lobarvalues(1)=0;
lobarvalues(2)=0;
lobarvalues(3)=0;
lobarvalues(4)=0;
lobarvalues(5)=0;
dicecoeff=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load the CT, Q and V images 
[LRCT_image,LRCTsize,dicomsize]=dicomload(imagepathCT);
[Q_image,SPECTsize, SPECTdicomsize]=dicomload(imagepathQ);
[V_image,SPECTsize,SPECTdicomsize]=dicomload(imagepathV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract lungs

% tidy up image
se = strel('disk', 20); 
Iobr = imerode(LRCT_image, se);
Iobr = imreconstruct(Iobr, LRCT_image);
Iobrd = imdilate(Iobr, se);
Iobrd = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
clear Iobr;

% invert image 
Iobrd = imcomplement(Iobrd);

% apply minimum threshold to image
Iobrd(Iobrd<lungmaxthreshold)=0;
if displaydebugimages >4
    figure('Name', 'Thresholded');
    imshow3Dfull(Iobrd);
    keyboard;
end;

% find optimal threshold using Otsu's method
values = [0 1];
thresh=multithresh(Iobrd);  
Iobrd = imquantize(Iobrd, thresh, values);                                              

% set all values on first slice where CT appears = 1
Iobrd(:,:,startsliceCT) = 1;
if displaydebugimages >4
    figure('Name', 'Binarised');
    imshow3Dfull(Iobrd);
    keyboard;    
end;

% dilate and fill in 2D then 3D
se = strel('sphere', 2);
Iobrd = imdilate(Iobrd, se);
Iobrcbrfill=Iobrd;
for i=1:LRCTsize(3)
    Iobrcbrfill(:,:,i) = imfill(Iobrcbrfill(:,:,i), 'holes');
end;
if displaydebugimages >4
    figure('Name', 'Filled');
    imshow3Dfull(Iobrcbrfill);
    keyboard;    
end;
Iobrcbrfill=imfill(Iobrcbrfill, 'holes');
Iobrcbrfill(:,:,startsliceCT) = 0;

% obtain lung by subtracting region from dilated and filled region then
% dilate and fill
lungROI=Iobrcbrfill-Iobrd;
clear Iobrcbrfill;
clear Iobrd;
lungROI = imdilate(lungROI, se);
for i=1:LRCTsize(3)
    lungROI(:,:,i) = imfill(lungROI(:,:,i), 'holes');
end;
lungROI=imfill(lungROI, 'holes');
Iobrcbrfill(:,:,startsliceCT) = 0;
if displaydebugimages >4
    figure('Name', 'Lung ROI');
    imshow3Dfull(lungROI);
    keyboard;    
end;

% erode lung
for i = 1:numerosions
    lungROI = imerode(lungROI, se);
end
if displaydebugimages >4
    figure('Name', 'Eroded Lung ROI');
    imshow3Dfull(lungROI);
    keyboard;    
end;

%%% Separate lungs using connected components
CC = bwconncomp(lungROI, conn);
numPixels = cellfun(@numel,CC.PixelIdxList);
s = regionprops(CC,'Centroid');
centroids = cat(1, s.Centroid);
top1=10000;
top2=10000;
numcomp=size(numPixels,2);
numbigcomp=0;
for i=1:numcomp
    if (numPixels(i)>mincomponentsize)
        numbigcomp=numbigcomp+1; 
        bigcentroids(numbigcomp,1)=centroids(i, 1);
        bigcentroids(numbigcomp,2)=centroids(i, 2);
        bigcentroids(numbigcomp,3)=centroids(i, 3);
        bigcentroids(numbigcomp,4)=i;
    end;
end;

% find 2 large conn comp with min y value
[big1,big1I]=min(bigcentroids(:,2));
bigcentroids(big1I,2)=500;
[big2,big2I]=min(bigcentroids(:,2));

% determine which is right and which is left
leftlungROI=lungROI;
rightlungROI=lungROI;
leftlungROI(:,:,:)=0;
rightlungROI(:,:,:)=0;
if bigcentroids(big1I,1)<bigcentroids(big2I,1)
    % big2 is left lung
    leftlungROI(CC.PixelIdxList{bigcentroids(big2I,4)})=1;
    rightlungROI(CC.PixelIdxList{bigcentroids(big1I,4)})=1;
else
    % big2 is right lung
    leftlungROI(CC.PixelIdxList{bigcentroids(big1I,4)})=1;
    rightlungROI(CC.PixelIdxList{bigcentroids(big2I,4)})=1;    
end;

for i = 1:numerosions
    leftlungROI = imdilate(leftlungROI, se);
    rightlungROI = imdilate(rightlungROI, se);
end

% create image of both lungs for checking
lungROI=leftlungROI+rightlungROI+rightlungROI;
if displaydebugimages >3
    figure('Name', 'Separated lungs');
    imshow3Dfull(lungROI);
    keyboard;
end;

writedicom( outputlabels, imagepathCT, lungROI, dicomsize, 'Matlab processed by EL', 'Separated lungs');

%%%  reslice LRCT to SPECT space
ny= SPECTsize(2);
nx= SPECTsize(1);
nz= SPECTsize(3);
[y x z]=ndgrid(linspace(1,size(leftlungROI,1),ny), linspace(1,size(leftlungROI,2),nx), linspace(1,size(leftlungROI,3),nz));
leftlungROI=interp3(leftlungROI, x,y,z, 'cubic');
rightlungROI=interp3(rightlungROI, x,y,z, 'cubic');
lungROI=interp3(lungROI, x,y,z, 'cubic');

if displaydebugimages >2
    figure('Name', 'L Lung ROI (before reg)');
    imshow3Dfull(leftlungROI);
    figure('Name', 'R Lung ROI (before reg)');
    imshow3Dfull(rightlungROI);
end;
    
[templatelobe_1,templatelobe_2,templatelobe_3,templatelobe_4,templatelobe_5,templatelungL,templatelungR,lungatlas ]= ...
   createLobes(template);

leftoverlapbefore=calcdicecoeff(leftlungROI, templatelungL);
rightoverlapbefore=calcdicecoeff(rightlungROI, templatelungR);

[leftlungregistered, lobe_1,lobe_2,lobe_3]=registerLungLobes(templatelungL, leftlungROI, SPECTsize, useSDFimages, ...
    lefttranslateflag, leftrigidflag, leftsimilarityflag, leftaffineflag, leftnonlinearflag, ...
    templatelobe_1,templatelobe_2,templatelobe_3);

%[rightlungregistered, lobe_4,lobe_5]=registerLungLobes(templatelungR, rightlungROI, SPECTsize, useSDFimages, ...
%    righttranslateflag, rightrigidflag, rightsimilarityflag, rightaffineflag, rightnonlinearflag, ...
%    templatelobe_4,templatelobe_5);

leftdicecoeff=calcdicecoeff(leftlungROI, leftlungregistered);
rightdicecoeff=calcdicecoeff(rightlungROI, rightlungregistered);

lobarvalues=measureLobarValues(lobe_1, lobe_2, lobe_3, lobe_4, lobe_5, V_image, Q_image);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Perform registration      %%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = registerLungLobes(varargin)

templatelung = varargin{1};
lungROI = varargin{2};
SPECTsize = varargin{3}; 
useSDFimages = varargin{4}; 
translate = varargin{5};
rigid = varargin{6};
similarity = varargin{7};
affine = varargin{8};
nonlinear = varargin{9};

global displaydebugimages;
resize=0;
numLobes= nargin-9;
pyramid1its=1000;
pyramid2its=1000;
pyramid3its=50;

registered_lobes=templatelung;
for i = 1:numLobes
    varargout{i+1} = varargin{i+9};
end;

if useSDFimages == 1    
    % convert images to SDF representation    
    templatelungSDF=convertBinToSDF(templatelung);
    lungROISDF=convertBinToSDF(lungROI);
    moving_image = templatelungSDF;
    fixed_image = lungROISDF;
else
    moving_image = templatelung;
    fixed_image = lungROI;
end;     
    
if displaydebugimages >2
    figure('Name', 'Moving image');
    imshow3Dfull(moving_image);
    figure('Name', 'Fixed image');
    imshow3Dfull(fixed_image);
end;

if displaydebugimages >0
    displayFusedImagePair(moving_image, fixed_image, SPECTsize, 'Fused pair before');
end;

if (nonlinear == 1)
    [D,registered_lobes]= imregdemons(moving_image, fixed_image, [pyramid1its,pyramid2its,pyramid3its], 'AccumulatedFieldSmoothing', 2);
    if displaydebugimages >2
        figure('Name', 'Moving image');
        imshow3Dfull(moving_image);
        figure('Name', 'Register lobes');
        imshow3Dfull(registered_lobes);
    end;
    if displaydebugimages >0
        displayFusedImagePair(registered_lobes, fixed_image, SPECTsize, 'Fused pair after');
    end;
    for i=1:numLobes   
        varargout{i+1} = imwarp(varargin{i+9}, D);
    end;
elseif (translate + rigid + similarity + affine > 0)    
    % Set optimizer and metric settings
    [optimizer,metric]= imregconfig('monomodal');

    % set that works for patient 4
    %optimizer.GradientMagnitudeTolerance = 1.000000e-09
    %optimizer.MinimumStepLength = 1.000000e-07
    %optimizer.MaximumStepLength = 1e-4
    %optimizer.MaximumIterations = 600
    %optimizer.RelaxationFactor = 0.99

    optimizer.GradientMagnitudeTolerance = 1.000000e-08
    optimizer.MinimumStepLength = 1.000000e-10
    optimizer.MaximumStepLength = 1e-7
    optimizer.MaximumIterations = 5000
    optimizer.RelaxationFactor = 0.999

    RA = imref3d(size(fixed_image));
    initform = affine3d();

    if translate==1
        % Perform 1st step of registration - Transformation
        tform1 = imregtform(moving_image,fixed_image,'translation',optimizer,metric, 'DisplayOptimization', true);
        temp1registered_lobes = imwarp(moving_image, tform1, 'OutputView', RA);
        temp1registered_lobes(temp1registered_lobes>0)=1;
        if displaydebugimages >0
            displayFusedImagePair(temp1registered_lobes, fixed_image, SPECTsize);
        end;
        overlapreg1=calcdicecoeff(lungROI, temp1registered_lobes, 'Fused pair after'); 
        dicecoeff=overlapreg1;
        initform=tform1;
    end;

    if rigid==1
        % Perform step of registration - 'rigid'
        tform2 = imregtform(moving_image,fixed_image,'rigid',optimizer,metric, 'InitialTransformation', initform, 'DisplayOptimization', true, 'PyramidLevels',4);
        temp2registered_lobes = imwarp(moving_image, tform2, 'OutputView', RA);    
        temp2registered_lobes(temp2registered_lobes>0.5)=1;
        if displaydebugimages >0
            displayFusedImagePair(temp2registered_lobes, fixed_image, SPECTsize, 'Fused pair after');
        end;
        overlapreg2=calcdicecoeff(lungROI, temp2registered_lobes);
        dicecoeff=overlapreg2;
        initform=tform2;
    end;
   
    if similarity==1 
        % Perform step of registration - 'similarity'
        tform3 = imregtform(moving_image,fixed_image,'affine',optimizer,metric, 'DisplayOptimization', true, 'PyramidLevels',4);
        temp3registered_lobes = imwarp(moving_image, tform3,  'OutputView', RA);
        if displaydebugimages >1
            figure('Name', 'Moving before similarity reg');
            imshow3Dfull(moving_image);
            figure('Name', 'fixed before similarity reg');
            imshow3Dfull(fixed_image);
            figure('Name', 'Moving after similarity reg');
            imshow3Dfull(temp3registered_lobes);
        end;
        temp3registered_lobes(temp3registered_lobes>0)=1;
        if displaydebugimages >0
            displayFusedImagePair(temp3registered_lobes, fixed_image, SPECTsize, 'Fused pair after');
        end;
        overlapreg3=calcdicecoeff(lungROI, temp3registered_lobes);
        dicecoeff=overlapreg3;
        initform=tform3;
        
    end;

    if affine==1
        % Perform step of registration - 'affine'
        if resize==1
            moving_image=resize(moving_image, size(moving_image).*[1/2,1/2,1/2]);
            fixed_image=resize(fixed_image, size(fixed_image).*[1/2,1/2,1/2]);
        end;
        optimizer.MaximumIterations = 100;
        tform = imregtform(moving_image,fixed_image,'affine',optimizer,metric, 'InitialTransformation', initform, 'DisplayOptimization', true, 'PyramidLevels',4);
        RA = imref3d(size(moving_image));
        registered_lobes = imwarp(moving_image, tform, 'OutputView', RA);
        if resize==1
            fixed_image = resize(fixed_image, size(fixed_image).*[2,2,2]);
            registered_lobes = resize(registered_lobes, size(moving_image).*[2,2,2]);
        end;
        
        if displaydebugimages >1
            figure('Name', 'Moving before similarity reg');
            imshow3Dfull(moving_image);
            figure('Name', 'fixed before similarity reg');
            imshow3Dfull(fixed_image);
            figure('Name', 'Moving after affine reg');
            imshow3Dfull(registered_lobes);
            figure('Name', 'Fixed after affine reg');
            imshow3Dfull(fixed_image);
            keyboard;
        end;
        if displaydebugimages > 0
            displayFusedImagePair(registered_lobes, fixed_image, SPECTsize,'Fused pair after');
        end;
        overlapreg4=calcdicecoeff(fixed_image, registered_lobes);
        dicecoeff=overlapreg4;
    end;
    for i=1:numLobes   
        varargout{i+1} = imwarp(varargin{i+9}, tform, 'OutputView', RA);
         if displaydebugimages > 0
            displayFusedImagePair(registered_lobes, fixed_image, SPECTsize, 'Fused pair after');
        end;
    end;
end;

if useSDFimages == 1
    registered_lobes(registered_lobes>0)=1;
    dicecoeff=calcdicecoeff(lungROI, registered_lobes);
    if displaydebugimages >1
            figure('Name', 'Moving after reg (binarised)');
            imshow3Dfull(registered_lobes);
            figure('Name', 'Fixed after reg (binarised)');
            imshow3Dfull(lungROI);
    end;
end;
varargout{1}=registered_lobes;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Measure Lobar Values      %%%%%%%%%%%%%%%%%%%%%%%%%
function [lobarvalues] = measureLobarValues(lobe_1, lobe_2, lobe_3, lobe_4, lobe_5, V_image, Q_image)

% % Lobar SPECT data
V_lobe1 = lobe_1.*V_image;
V_lobe2 = lobe_2.*V_image;
V_lobe3 = lobe_3.*V_image;
V_lobe4 = lobe_4.*V_image;
V_lobe5 = lobe_5.*V_image;
Q_lobe1 = lobe_1.*Q_image;
Q_lobe2 = lobe_2.*Q_image;
Q_lobe3 = lobe_3.*Q_image;
Q_lobe4 = lobe_4.*Q_image;
Q_lobe5 = lobe_5.*Q_image;

% Sum of lobe values is greater than one because they overlap sometimes
lobe_value1 = sum(sum(sum(V_lobe1))) ;
lobe_value2 = sum(sum(sum(V_lobe2))) ;
lobe_value3 = sum(sum(sum(V_lobe3))) ;
lobe_value4 = sum(sum(sum(V_lobe4))) ;
lobe_value5 = sum(sum(sum(V_lobe5))) ;

total_emission = lobe_value1 + lobe_value2 + lobe_value3 + lobe_value4 + lobe_value5;
 
lobarvalues(1)=lobe_value1 / total_emission;
lobarvalues(2)=lobe_value2 / total_emission;
lobarvalues(3)=lobe_value3 / total_emission;
lobarvalues(4)=lobe_value4 / total_emission;
lobarvalues(5)=lobe_value5 / total_emission;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write dicom images
function [] = writedicom( outputimagepath, templateimagepath, img1, size, seriesdesc, imageID)

img1=cast(img1,'uint16');
metadata = dicominfo(templateimagepath);
img1=reshape(img1,size);

uid = dicomuid;
metadata.SeriesInstanceUID = uid;
metadata.SeriesDescription=seriesdesc;

metadata.ImageID = imageID;
dicomwrite(img1, outputimagepath, metadata,'CreateMode', 'copy'); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Measure Dice coefficients %%%%%%%%%%%%%%%%%%%%%%
function [dicecoeff] = calcdicecoeff( img1, img2)

%%%%%% Lung ROI Dice coefficient
intersection=img1&img2;
volimg1=sum(img1(:));
volimg2=sum(img2(:));
volintersect=sum(intersection(:));
dicecoeff=2*volintersect/(volimg1+volimg2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Convert binary image to SDF %%%%%%%%%%%%%%%%%%%%%%
function [sdfImage] = convertBinToSDF( binImage )

 binImage=binImage>0;
 inSDF=bwdist(binImage);
 comp=imcomplement(binImage);
 outSDF=bwdist(comp);
 sdfImage=outSDF-inSDF;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Load dicom images %%%%%%%%%%%%%%%%%%%%%%%%
function [image, imagesize, dicomsize]=dicomload(imagepath)

global displaydebugimages;
info = dicominfo(imagepath);
image = dicomread(info);
dicomsize=size(image);
image = squeeze(image);
imagesize=size(image);
if displaydebugimages > 5
    figure('Name', imagepath);
    imshow3Dfull(LRCT_image);
    keyboard;
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Display fused image pair %%%%%%%%%%%%%%%%%%%%%%%%
function displayFusedImagePair(image1, image2, size, name)

global displaydebugimages;
image1_col=image1(:);
image2_col=image2(:);
fused_image = imfuse(image1_col, image2_col,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
fused_image = reshape(fused_image,size(1),size(2),size(3), 3);
SliceBrowser(fused_image);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      Produce lobar map       %%%%%%%%%%%%%%%%%%%%%%%
function [templatelobe_1,templatelobe_2,templatelobe_3,templatelobe_4,templatelobe_5,templatelungL,templatelungR,lungatlas ]= ...
   createLobes(template)

global displaydebugimages;
displaytemplatelobeimages=0;

%%% Load the template image
fid = fopen(template,'r');
lobes = fread(fid,5* 128* 128* 128,'*int8');
fclose(fid);
lobes = reshape(lobes, 5,128,128,128);
if displaytemplatelobeimages == 1
    figure('Name', 'Lobe 1');
    imshow3Dfull(squeeze(lobes(1,:,:,:)));
    figure('Name', 'Lobe 2');
    imshow3Dfull(squeeze(lobes(2,:,:,:)));
    figure('Name', 'Lobe 3');
    imshow3Dfull(squeeze(lobes(3,:,:,:)));
    figure('Name', 'Lobe 4');
    imshow3Dfull(squeeze(lobes(4,:,:,:)));
% orientation of image is lobes(lobenumber, y, x, slice)
    figure('Name', 'Lobe 5');
    imshow3Dfull(squeeze(lobes(5,:,:,:)));
    keyboard;    
end;   
   
templatelobe_1 = uint16(squeeze(lobes(1,:,:,:)));
templatelobe_2 = uint16(squeeze(lobes(2,:,:,:)));
templatelobe_3 = uint16(squeeze(lobes(3,:,:,:)));
templatelobe_4 = uint16(squeeze(lobes(4,:,:,:)));
templatelobe_5 = uint16(squeeze(lobes(5,:,:,:)));

% Map of right template lung for registration purposes
templatelungR = templatelobe_1+templatelobe_2+templatelobe_3;
templatelungR(templatelungR>0)=1;

% Map of left template lung for registration purposes
templatelungL = templatelobe_4+templatelobe_5;
templatelungL(templatelungL>0)=1;

% Map of template lungs (L+R) for registration purposes
lungatlas = templatelungL + templatelungR;

if displaydebugimages >2
    figure('Name', 'L Lung atlas (before reg)');
    imshow3Dfull(templatelungL);
    figure('Name', 'R Lung atlas (before reg)');
    imshow3Dfull(templatelungR);
    keyboard;
end;

end

