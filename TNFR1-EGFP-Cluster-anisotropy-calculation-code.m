% This code calculates fluorescence anisotropy for small objects of micron size based on images.
% We use this code to calculate Fl anisotropy of membrane-associated TNFR1-EGFP clusters.
% Imaging was performed with a Zeiss AxioObserver Z1 epifluorescence microscope equipped 
% with a 63x oil immersion objective (NA 1.4) and a mercury arc lamp (HXP 120 V) as the 
% light source. Horizontally polarized excitation light was achieved by passing the light 
% through a linear polarizer (ThorLabs) in the excitation path. In the emission path, 
% a polarizing beam splitter (DV2, Photometrics) was utilized to divide the resulting 
% polarized fluorescence signal into parallel and perpendicular polarizations. 
% The fluorescence from these two polarizations was simultaneously captured by a CMOS camera
% (Hamamatsu Orca Flash 4.0 C13440), with one half of the image representing the parallel channel (I_∥)
% and the other half representing the perpendicular channel (I_⊥). 
% Due to deformations in the optical path and thermal fluctuations, there was partial overlap between
% the two halves (I_∥, I_⊥). To correct this, image registration of each cluster of TNFR1-EGFP in the
% two channels was performed by intensity cross-correlation using the particle image velocimetry (PIV)
% method followed by MATLAB's mono-modal intensity-based registration tool. 
% Subsequently, the fluorescence anisotropy (r) in each pixel of the clusters was calculated 
% using the equation r = (I_∥- g I_⊥)/(I_∥+ 2 g I_⊥ ).Where g is the instrumental correction factor 
% (G-factor, I_∥ /I_⊥ ), determined by measuring 100 nM fluorescein in water in each pixel.  
% Here left side image is perperdicular ch(FIXED) and Right side image is Parallel ch(MOVING)
% Before running the program please check whether the matlab function 
% registerImage.m and Gfactor image and this code are there in same directory. 
% Check the x-y shift of cluster for an image. Upon running the code first you
% have to chose the file and then take background(BG) from left side image, where
% there is no cell image, then take the cell of interest and the cluster you want.
% Image will be saved as Tif in working directory.

close all

PixShift = 3; % for PIV
c0 = 1; %cell no
% from raw image try to evaluate how much shift in between 
%parrallel and perpendicluar channel is there in x and y for a certain cluster
dx = 1023; dy =3;  % x-shift and y shift ;change here

%input file
[file,path] = uigetfile('*.tif*');
f1 = fullfile(path,file);
im_info = imfinfo(f1);
Nimg = size(im_info,1);
Gfactor = imread('Gfactor_63x_Hpol_polarizer.tif');
for n_img =1:Nimg % for z stack image start n_img from focus basal plane of your interest
    a = imread(f1,n_img);
    % Background correction 
    [left_back,rect0] = imcrop(a); %crop image from no cell region
    back_left =mean(mean(left_back)); %background correction
   % back_left = 0;
    right_back = imcrop(a,[rect0(1)+ dx rect0(2)+dy rect0(3) rect0(4)]);
    back_right = mean(mean(right_back));%background correction
   % back_right =0;
   
    %cluster formation 
    b=imadjust(a,[0.05 0.50]);
    [Left1,rect] = imcrop(b);
    Left = imcrop(a,[rect(1) rect(2) rect(3) rect(4)]);
    Right = imcrop(a,[rect(1)+dx rect(2)+dy rect(3) rect(4)]);
    Gfact_Left = imcrop(Gfactor,[rect(1) rect(2) rect(3) rect(4)]);
    Gfact_Right = imcrop(Gfactor,[rect(1)+ dx rect(2)+ dy rect(3) rect(4)]);
    imshow(Left);
    

    it = input('Initial start value of cluster : staring roi = ');
    t = input('How many focused clusters you want to calculate at a time ? : no of roi = ');
    for roi = it:t
        rect_sub = cell(1,t-it);
        [Left_sub_cluster1,rect_sub{roi}] = imcrop(Left1);
        Left_sub_cluster = imcrop(Left,rect_sub{1,roi});
        Right_sub_cluster = imcrop(Right, rect_sub{1,roi});
        Gfact_sub_Left = imcrop(Gfact_Left,rect_sub{1,roi});
        Gfact_sub_Right = imcrop(Gfact_Right,rect_sub{1,roi});
        

% PIV
        piv=[];
        PIV_intensity =[];
        
        for i = -PixShift : PixShift 
            for j = -PixShift : PixShift
                x1 = rect_sub{1,roi}(1)+i;
                y1 = rect_sub{1,roi}(2)+j;
                move_Right_sub = imcrop(Right,[x1 y1 rect_sub{1,roi}(3) rect_sub{1,roi}(4)]);
                
                multiply_left_right = immultiply(double(Left_sub_cluster),double(move_Right_sub));
                sum_intensity = sum(sum(multiply_left_right));
                piv =[x1,y1,sum_intensity];
                PIV_intensity = [PIV_intensity; piv];
              
           
            end
        end
        [q,w]= size(PIV_intensity);
        maxValue1 = max(PIV_intensity(:)); 
        [rM ,cM] = find(PIV_intensity == maxValue1);
        xMaxShift = PIV_intensity(rM, cM-2); yMaxShift = PIV_intensity(rM, cM-1);
        registered_Right_sub_cluster = imcrop(Right,[xMaxShift yMaxShift rect_sub{1,roi}(3) rect_sub{1,roi}(4)]);
        r_rsc = double(registered_Right_sub_cluster );
        r_lsc = double(Left_sub_cluster);

%%%%% MATLAB REGISTRATION :Intensity Based Monomodal Registration
        MOVING = double(r_rsc);
        FIXED = double(r_lsc);
    
        registerImages(MOVING,FIXED);
        reg_img = ans.RegisteredImage;
   
   
%%%Anisotropy Calculation with MATLAB Registered image
   
        Parallelx= double(reg_img) ;
        [nx,ny] = size(Parallelx);
        Parallel = imsubtract(double(Parallelx),double(back_right)); %background correction
        Perpendicularx= double(r_lsc);
        Perpendicular =  imsubtract(double(Perpendicularx),double(back_left)); %background correction
    
        Gfactor_m = imdivide(double(Gfact_sub_Right),double(Gfact_sub_Left));
        Perp_Gfact = immultiply((Perpendicular), double(Gfactor_m));  % G-factor multiplication 
        TotalIntensity_perp = immultiply(double(Perp_Gfact),2);
        TotalIntensity_mat = imadd(double(Parallel),double(TotalIntensity_perp));
        Polarisation = imsubtract(double(Parallel),double(Perp_Gfact));
        Anisotropy_mat =  imdivide(double(Polarisation),double(TotalIntensity_mat));
 
        % name of files of anisotropy and intensity image
        ani_tif = strcat('Ani_ctrl-HeLa-TNFR1GFP_cel0',num2str(c0),'_Clus',num2str(roi),'_xsh',num2str(dx),'_ysh',num2str(dy),'.tif');
        totInt_tif = strcat('TotInt_ctrl-HeLa-TNFR1GFP_cel0',num2str(c0),'_Clus',num2str(roi),'_xsh',num2str(dx),'_ysh',num2str(dy),'.tif');
       
        %save anisotropy Image as a 32bit .tif using TIFF class
 
            a_peak = Tiff(char(ani_tif),'w');

            setTag(a_peak,'Photometric',Tiff.Photometric.MinIsBlack);
            setTag(a_peak,'Compression',Tiff.Compression.None);
            setTag(a_peak,'BitsPerSample',32);
            setTag(a_peak,'SamplesPerPixel',1);
            setTag(a_peak,'SampleFormat',Tiff.SampleFormat.IEEEFP);
            setTag(a_peak,'ImageLength',nx);
            setTag(a_peak,'ImageWidth',ny);
            setTag(a_peak,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    
            write(a_peak,single(Anisotropy_mat));
            a_peak.close();
        
            
          %save TotalInt Image as a 32bit .tif using TIFF class   
            t_peak = Tiff(char(totInt_tif),'w');

            setTag(t_peak,'Photometric',Tiff.Photometric.MinIsBlack);
            setTag(t_peak,'Compression',Tiff.Compression.None);
            setTag(t_peak,'BitsPerSample',32);
            setTag(t_peak,'SamplesPerPixel',1);
            setTag(t_peak,'SampleFormat',Tiff.SampleFormat.IEEEFP);
            setTag(t_peak,'ImageLength',nx);
            setTag(t_peak,'ImageWidth',ny);
            setTag(t_peak,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    
            write(t_peak,single(TotalIntensity_mat)); 
            t_peak.close();
           
    end
   
end
