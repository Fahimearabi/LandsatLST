%Estimation of land surface temperature by Landsat 9
%This code calculates the land surface temperature using the Landsat 9 satellite image.

%%!!Code execution steps
% Bands 4, 5, 10 and 11 of the Landsat 8 satellite are placed in one folder.
% Name them only as numbers and in tiff format.
% In this folder, the output file, which is the image of the surface temperature, is saved.

% This code was written by Fahime Arabi. if you have any questions about it, I will answer you with the following email:
%Fahimearabi1993@gmail.com
clc
clear all
close all
% -----------------------load images landsat8---------------------------
[B4,RefMatrx]=geotiffread('B4.TIF');
InfoB4=geotiffinfo('B4.tif');
[B5,RefMatrx]=geotiffread('B5.tif');
[B10,RefMatrx]=geotiffread('B10.tif');
[B11,RefMatrx]=geotiffread('B11.tif');
% --------------------------------------------------
B4=double(B4);
B5=double(B5);
B10=double(B10);
B11=double(B11);
% -----------------------Radiance---------------------------
radianceb10=(3.8000E-04*B10)+0.1;
radianceb11=(3.4900E-04*B11)+0.1;
% -----------------------CONSTANT---------------------------
k2b10=1329.2405;
k2b11=1198.3494;
k1b10=799.0284;
k1b11=475.6581;
% -----------------------Brightness_Temperature---------------------------
btb10=k2b10./log((k1b10./double(radianceb10))+1);
btb11=k2b11./log((k1b11./double(radianceb11))+1);
% -----------------------NDVI---------------------------
ndvi=(B5-B4)./(B5+B4);
[x,y,d]=size(ndvi);
s(x,y)=0;

for i=1:x;
    for j=1:y;
        
       if ndvi(i,j)>=0 & ndvi(i,j)<0.05;
           s(i,j)=ndvi(i,j);
       else
           s(i,j)=0;
       end
    end
end
b=nonzeros(s);
ndvis=mean(b);


v(x,y)=0;
for i=1:x;
    for j=1:y;
        
       if ndvi(i,j)>=0.2;
           v(i,j)=ndvi(i,j);
       else
           v(i,j)=0;
       end
    end
end

b1=nonzeros(v);
ndviv=mean(b1);
% --------------------------------------------------
fvc=(ndvi-ndvis)/(ndviv-ndvis);
% --------------------------------------------------
esb10=0.971;
esb11=0.977;
evb10=0.987;
evb11=0.989;
% -----------------------Emissivity---------------------------
e10=esb10*(1-fvc)+evb10*fvc;
e11=esb11*(1-fvc)+evb11*fvc;
e=(e10+e11)/2;
deltae=e10-e11;


% -----------------------CONSTANT---------------------------
w=0.013;
c0=0.286;
c1=1.376;
c2=0.183;
c3=54.3;
c4=2.238;
c5=129.2;
c6=16.4;
% -----------------------LSTSW1-SKOKOVIC--------------------------
lst=btb10+c1*(btb10-btb11)+c2*(btb10-btb11).^2+double((c3+c4*w)*(1-e)+(c5+c6*w)*deltae+c0);
lstSW1=lst-273.15;
disp('LSTSW1 executed successfully');
geotiffwrite('lstSW1',lstSW1,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);