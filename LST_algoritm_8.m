%Estimation of land surface temperature by Landsat 8
%This code calculates the land surface temperature using the Landsat 8 satellite image.
% In this code, the land surface temperature is calculated using single channel and separate window methods.

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
[B4,RefMatrx]=geotiffread('4.tif');
InfoB4=geotiffinfo('4.tif');
[B5,RefMatrx]=geotiffread('5.tif');
[B10,RefMatrx]=geotiffread('10.tif');
[B11,RefMatrx]=geotiffread('11.tif');
B4=double(B4);
B5=double(B5);
B10=double(B10);
B11=double(B11);
%% -----------------------Radiance---------------------------
radianceb10=(0.0003342*B10)+0.1;
radianceb11=(0.0003342*B11)+0.1;
%% -----------------------CONSTANT---------------------------
k2b10=1321.0789;
k2b11=1201.1442;
k1b10=774.8853;
k1b11=480.8883;
%% -----------------------Brightness_Temperature---------------------------
btb10=k2b10./log((k1b10./double(radianceb10))+1);
btb11=k2b11./log((k1b11./double(radianceb11))+1);
%% -----------------------NDVI---------------------------
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
%% --------------------------------------------------
fvc=(ndvi-ndvis)/(ndviv-ndvis);
%% --------------------------------------------------
esb10=0.971;
esb11=0.977;
evb10=0.987;
evb11=0.989;
% -----------------------Emissivity---------------------------
e10=esb10*(1-fvc)+evb10*fvc;
e11=esb11*(1-fvc)+evb11*fvc;
e=(e10+e11)/2;
deltae=e10-e11;


%%-----------------------CONSTANT---------------------------
w=0.013;
c0=0.286;
c1=1.376;
c2=0.183;
c3=54.3;
c4=2.238;
c5=129.2;
c6=16.4;

%%  -----------------------LSTSW1-SKOKOVIC--------------------------
lst=btb10+c1*(btb10-btb11)+c2*(btb10-btb11).^2+double((c3+c4*w)*(1-e)+(c5+c6*w)*deltae+c0);
lstSW1=lst-273.15;
disp('LSTSW1 executed successfully');

geotiffwrite('lstSW1',lstSW1,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%% -------------LSTSW2-KERR-----------
    MIN_NDVI=min(min(ndvi));
    Max_NDVI=max(max(ndvi));
    pv=(ndvi-(MIN_NDVI))./((Max_NDVI)-(MIN_NDVI));
    LSTSW2=(btb10.*(0.5*pv+3.1)+btb11.*(-0.5*pv-2.1)-5.5*pv+3.1)-273.15;
    disp('LSTSW2 executed successfully'); 
 geotiffwrite('LSTSW2',LSTSW2,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
   
 %--------------LSTSW3-MCCLAIN-------------------
 LSTSW3= (1.035* btb10+3.046*(btb10-btb11)-10.93)-273.15;
disp('LSTSW3 executed successfully');
geotiffwrite('LSTSW3',LSTSW3,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag); 
%--------------LSTSW4-PRATA-------------------
t0=310;
LSTSW4=3.45*((btb10-t0)./e10)-2.54*((btb11-t0)./e11)+40*((1-e10)./e10)+t0-273.15;
disp('LSTSW4 executed successfully');
geotiffwrite('LSTSW4',LSTSW4,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
%--------------------LSTSW5-PRICE----------------------
LSTSW5=((btb10+3.33*(btb10-btb11)).*((5.5-e10)./4.5)+0.75*btb11.*deltae)-273.15;
disp('LSTSW5 executed successfully');
geotiffwrite('LSTSW5',LSTSW5,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%-------------------LSTSW6-SOBRINO-------------------
LSTSW6=btb10+1.06*(btb10-btb11)+0.46*((btb10-btb11).^2)+53*(1-e10)-53*(e10-e11)-273.15;
disp('LSTSW6 executed successfully');
geotiffwrite('LSTSW6',LSTSW6,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%------------------LSTSW7-ULIVIERII-------------
LSTSW7=(btb10+3*(btb10-btb11)-52.45*e+51.57)-273.15;
disp('LSTSW7 executed successfully');
geotiffwrite('LSTSW7',LSTSW7,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%-------------LSTSW8-VIDAL-----------------------
LSTSW8=(3.78.*btb10)-(2.78.*btb11)+(50.*(1-(e10./e11)))-(300.*(deltae./e10))-273.15;
disp('LSTSW8 executed successfully');
geotiffwrite('LSTSW8',LSTSW8,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
%---------------LSTSW9-COLL----------------------------
LSTSW9=(0.39*(btb10.^2)+2.3*btb10-0.78.*btb10.*btb11-1.34.*btb11+0.39*(btb11.^2)+0.56)-273.15;
disp('LSTSW9 executed successfully');
geotiffwrite('LSTSW9',LSTSW9,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%--------------LSTSWCR10-------------------
Max_btb11=max(max(btb11));
Mean_btb11=mean(mean(btb11));
Mean_btb10=mean(mean(btb10));
R=((btb10-Mean_btb10).*(btb11-Mean_btb11))./((btb10-Mean_btb10).*(btb10-Mean_btb10));
tt=(e10./e11).*R;
w=(-13.412.*tt)+14.158;
LSTSWCR=(btb10+(1.376.*(btb11-btb11))+((0.183*(btb10-btb11)).^2)+((54.3+2.238.*w).*(1-e))+((129.2+16.4.*w).*deltae)+0.286)-273.15;
disp('LSTSWCR10 executed successfully')
geotiffwrite('LSTSWCR',LSTSWCR,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%----------------LSTSWCR11-feyzifrm------------------
w2=-9.674+0.653*tt+9.087.*(tt.*tt);
LSTSWCR2=(btb10+(1.376.*(btb11-btb11))+((0.183*(btb10-btb11)).^2)+((54.3+2.238.*w).*(1-e))+((129.2+16.4.*w2).*deltae)+0.286)-273.15;
disp('LSTLSTSWCR11 executed successfully')
geotiffwrite('LSTSWCR2',LSTSWCR2,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%----------------LSTSC1-ARTIS--------------
lne=(log10(e)/log10(exp(1)));
LSTSC1_10=(btb10./(1+10.8*(0.000001*(btb10./(1.438*0.01).*lne))))-273.15;
disp('LSTSc1_10 executed successfully');
LSTSC1_11=(btb11./(1+10.8*(0.000001*(btb11./(1.438*0.01).*lne))))-273.15;
disp('LSTSc1_11 01executed successfully');
geotiffwrite('LSTSC1_10',LSTSC1_10,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('LSTSC1_11',LSTSC1_11,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%% ------------LSTSC-Zhang-------------
LSTSC3_10=btb10./(((lne.*11.5.*btb10)./14380)+1)-273.15;
disp('LSTSc3_10 01executed successfully');
LSTSC3_11=btb11./(((lne.*11.5.*btb11)./14380)+1)-273.15;
disp('LSTSc3_11 01executed successfully');
geotiffwrite('LSTSC3_10',LSTSC3_10,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('LSTSC3_11',LSTSC3_11,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%% ---------------------LSTMW-QIN------------------------
C=e.*0.94;
D=(1-0.94).*(1+(1*e)*0.94);
a1=-62.7187;
a2=0.4339;
a3=315;
LSTMW1_10=((a1.*(1-C-D)+(a2.*(1-C-D)+C+D).*btb10-D.*a3)./(C))-273.15;
LSTMW1_11=((a1.*(1-C-D)+(a2.*(1-C-D)+C+D).*btb11-D.*a3)./(C))-273.15;
disp('LSTmw1_10 executed successfully')
disp('LSTmw1_11 executed successfully')
geotiffwrite('LSTMW1_10',LSTMW1_10,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('LSTMW1_11',LSTMW1_11,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

%% ----------------LSTMW2-RONGALL---------------------
LSTMW2_11=(btb11./1+(B11.*btb11./14380).*lne)-273.15;
LSTMW2_10=(btb10./1+(B10.*btb10./14380).*lne)-273.15;
disp('LSTmw2_10 executed successfully')
disp('LSTmw2_11 executed successfully')
geotiffwrite('LSTMW2_10',LSTMW2_10,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('LSTMW2_11',LSTMW2_11,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);

