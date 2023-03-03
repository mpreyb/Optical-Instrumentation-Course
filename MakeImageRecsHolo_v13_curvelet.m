% Useful Directories
clear all;close all;clc

%%This path must be changed as necessary
RADDIR='C:\Users\PcCristian\Downloads\PhaseMaps_DHI\';
CheminH=[RADDIR,'1_CCD_Acquisitions\plate2 3W 384x264 100000fps lenses15cm dr160 obj191cm size20cm 1825Hz 0.5V circle 100ratioPerfect exp1_996923 fiberless\'];
CheminR=[RADDIR,'1_CCD_Acquisitions\plate2_ref\'];
CheminA = [RADDIR,'2_Amplitudes\'];
CheminO = [RADDIR,'5_Phase_Diff\'];
CheminC = [RADDIR,'8_Phase_Diff_Filt_Curvelet\'];
CheminPDUW = [RADDIR,'7_Phase_Diff_Unwrap\'];
% mkdir(CheminA);
% mkdir(CheminPDF);
% mkdir(CheminPDUW);

imagefilesH = dir(CheminH);
imagefilesR = dir(CheminR);    

nfiles = length(imagefilesH);% Number of files found

disp(['Number of acquisitions: ',num2str(nfiles)]);

%Experimental parameters

NFFT=1024;
MFFT=NFFT;
px = 18.5;
py = px;
LAMBDA = 532;

%dr = input('Reconstruction distance mm: ');
dr = 160;

for m=3:nfiles
   
    tic
    %Loading acquisitions
    currentfilename = imagefilesH(m).name;
    H(:,:,m-2) = double(imread([CheminH,currentfilename]));
    
    Hm = H(:,:,m-2);
    
    if m == 3
        %figure(1);imagesc(Hm);colormap('gray');title('Hologram');
        %FFT_H = fftshift(fft2(fftshift(Hm)));
        %figure(1);imagesc(abs(log(FFT_H)));colormap('gray');title('FFT Hologram');
    end
   
    %%figure;imagesc(Hm);title([num2str(m-2),' hol']);axis off; axis equal
    %%colormap(gray);
    
    currentfilename = imagefilesR(m).name;
    R(:,:,m-2) = double(imread([CheminR,currentfilename]));
    
    Hm = R(:,:,m-2);
    
    %%figure;imagesc(Hm);title([num2str(m-2),' ref']);axis off; axis equal
    %%colormap(gray);
    
end

%
Rmean = mean(R,3);
%figure;imagesc(Rmean);title('Rmean');axis off; axis equal
%colormap(gray);.

% Amplitude calculation

for m=1:nfiles-2
%for m=72:72

    tic
    Hm = R(:,:,m);
     
    %figure;imagesc(Hm);title([num2str(m-2),' hologram']);axis off; axis equal
    %colormap(gray);
  
    Hm = H(:,:,m);
    Hp(:,:,m) = Hm - Rmean;
    
  
    %Numerical propagation   
    Complex(:,:,m) = TdeFresnel2(Hp(:,:,m),MFFT,NFFT,px,py,dr,LAMBDA,1e12,1,1);
    
     % Amplitude extraction
     Cm = Complex(:,:,m);
     A=abs(Cm)+eps; 
     %B = A;
     B=conv2(A,ones(20,20)/400,'same');
     %H_AMP=log10(B);
     H_AMP=abs(B);
     
     %figure;imagesc(H_AMP);title([num2str(m),' proper reconstructed hologram amplitude abs(Ar)']);axis off; axis equal
     %colormap(gray);
    
     if m == 28
        % Save this amplitude 
        MAXAmp = max(max(H_AMP));
        MINAmp = min(min(H_AMP));
        FileAmp = (H_AMP- MINAmp)/(MAXAmp-MINAmp);
        FileAmp = 255*FileAmp;
        imwrite(FileAmp, gray(256), [CheminA, int2str(dr),  'mm proper reconstructed hologram amplitude', '.tif']);
     end
     
   
    % Phase
    Phi(:,:,m) = angle(Complex(:,:,m));
    
    toc

end

 %% focal plane amplitude finding 
 %You don't need to run this part of the code is you know the actual
 %reconstruction distance
 
 for kk=100:2:200
         
     tic
     
     Hm = H(:,:,m);
     Hp(:,:,m) = Hm - Rmean;
     
     [Ar] = TdeFresnel2(Hp(:,:,m),MFFT,NFFT,px,py,kk,LAMBDA,1e12,1,1);
         
     A=abs(Ar)+eps;
     B=conv2(A,ones(10,10)/100,'same');
     H_AMP=log10(B);
     H_AMP=abs(H_AMP);
  
     figure;imagesc(abs(H_AMP));title([num2str(kk),' mm ']);axis off; axis equal
     colormap(gray);
     toc
     
     kk
 end

%% Phase calculation

for m = 1:nfiles-3
%for m = 1:3
%for m = 9:9

   
   % Phase difference
   Phi_diff(:,:,m) = Phi(:,:,m+1) - Phi(:,:,m);
   % Phase diff denoising
 
   
   PhiComplex(:,:,m) = exp(j.*Phi_diff(:,:,m));
   MAXPhase0= max(max(PhiComplex(:,:,m)));
   MINPhase0 = min(min(PhiComplex(:,:,m)));
   FilePhase0 = (PhiComplex(:,:,m)- MINPhase0)/(MAXPhase0-MINPhase0);
   FilePhase0 = 255*FilePhase0;
   imwrite(FilePhase0, gray(256), [CheminO, int2str(m),  ' PhaseDiff', '.tif']);
   
   %cvt filter by thresholding
   tic
   
   J = 4; 
   L = [3 4 4 5];
   
   % Noise estimation
   [wca wch wcv wcd] = dwt2( PhiComplex(:,:,m), 'db2' );
   sig = median(abs( wcd(:) ))/0.001;
   disp('noise estimation finished');
  
   cv = cvt(PhiComplex(:,:,m) , J, L, 1);
   disp('cvt finished');
   load cvt_th_4_3445_mean_lasl_zero.mat
   cth = cvt_llas_2_lasl( cvt_th );
   cth = cellmul( cth, 10*sig );
   nn = length( cth );
   cth{nn} = cellmul( cth{nn}, 10 );
   cv = cvt_llas_2_lasl( cv );
   yd = cellmul( cv, cellcompare( cellabs( cv ), cth ) );
   yd = cvt_lasl_2_llas( yd );
   disp('detection finished');
   PhiFiltered_curvelet(:,:,m) = icvt( yd, J, L, 1 );
   disp('icvt finished');
   tim = toc/60; disp('time = ');
   disp(tim);
   psnr_est = psnr( PhiComplex(:,:,m), PhiFiltered_curvelet(:,:,m) );
   
   
   
   %Save this phase difference filtered by cvt
   MAXPhase= max(max(PhiFiltered_curvelet(:,:,m)));
   MINPhase = min(min(PhiFiltered_curvelet(:,:,m)));
   FilePhase = (PhiFiltered_curvelet(:,:,m)- MINPhase)/(MAXPhase-MINPhase);
   FilePhase = 255*FilePhase;
   imwrite(FilePhase, gray(256), [CheminC, int2str(m),  ' PhaseDiffFiltered_curvelet', '.tif']);
   
   if m == 1
        figure;hold on;colormap(gray);imagesc(H_AMP);axis tight;

        %Determination of the Mask that will be used (circle)

        Npts=4;
        title(['Cliquer sur ',num2str(Npts),' bord du cercle du poinçon de l''objet']);
        hold on;
        for q1=1:Npts
            [XYG]=ginput(1);
            xg(q1)=XYG(1,1);yg(q1)=XYG(1,2);
            plot(xg(q1),yg(q1),'*r');
        end 

        % circle determination
        [Xcer,Ycer,Rcer]=CalculCercle(xg,yg);
        teta=linspace(0,2*pi,50);
        plot(Xcer+Rcer*cos(teta),Ycer+Rcer*sin(teta),'r',Xcer,Ycer,'or');

        Uo = Xcer;
        Vo = Ycer;
        Ru = Rcer;

        % Mask creation
        xe=[0:NFFT-1];% vecteur x
        ye=[0:MFFT-1];% vecteur y
        [U,V]=meshgrid(xe,ye);

        IF=find((U-Uo).^2+(V-Vo).^2 <= Ru^2);
        Mask=zeros(size(H_AMP));
        Mask(IF)=ones(size(IF)); 

        figure;colormap(gray);imagesc(Mask);axis tight;
    
        %unwrapping

        figure;hold on;colormap(gray);imagesc(PhiFiltered_curvelet(:,:,m).*Mask);axis tight;

        disp('Select 1 phase_known point on the wrapped phase map');
        [XYG]=ginput(1);
        XG=round(XYG(1,2));YG=round(XYG(1,1));     

   end


   toc
    
end

%% Unwrapping via CPULSI

for m = 1:nfiles-3
%for m = 21:2:41
%for m = 9:9

[phase_unwrap(:,:,m),phase_calibrate(:,:,m),N_iteration,t] = CPULSI(PhiFiltered_wft2f(:,:,m),Mask, 100,0.01,XG,YG,'true');

figure;colormap(gray);imagesc(phase_unwrap(:,:,m).*Mask);title([' phase unwrap ',num2str(m)]);axis tight;
%figure;colormap(gray);imagesc(phase_calibrate(:,:,m).*Mask);title([' phase_calibrate ',num2str(m)]); axis tight;

%Save this phase difference filtered by wft2f
MAXPhase= max(max(phase_unwrap(:,:,m).*Mask));
MINPhase = min(min(phase_unwrap(:,:,m).*Mask));
FilePhase = (phase_unwrap(:,:,m).*Mask- MINPhase)/(MAXPhase-MINPhase);
FilePhase = 255*FilePhase;
imwrite(FilePhase, gray(256), [CheminPDUW, int2str(m),  ' phase unwrap', '.tif']);

end

