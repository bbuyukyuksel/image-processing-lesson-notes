%%Matlab kullananlar için aşağıdaki satır yorum satırı olarak bırakılmalıdır.
%pkg load image;

clc; clear; close all;
%% İlk olarak 'Merhaba Dünya'
printf("Hello World \n\r");

%% Kullanıcıdan sisteme girdi almak için input
%% ve cout gibi yardıma yetişen disp

disp("Benim sayım 2,");
sayi = input('Bir sayi girsene : ');
disp(['Bak toplamları: ', num2str(sayi+ 2) , ' oldu.'])

%% Uykular için pause(time);
 
chess = zeros(10,10);
dolu = '#'
bos = '.'
for i=1:1:10
    for j=1:1:10
        if(mod(i,2) == 1) %Satırın tek oldugu degerlerde islem yap
            if(mod(j,2) == 0)
                %chess(i,j) = '#'
                chess(i,j) = dolu;
                printf("%c ", dolu);
            else
                chess(i,j) = bos;
                printf("%c ", bos)
            end
        else
            if(mod(j,2) == 1)
                chess(i,j) = dolu;
                printf("%c ", dolu);
            else
                chess(i,j) = bos;
                printf("%c ", bos);
            end
        end
    end
    printf("\n");
end

% Disp diyince sonucumuz neden böyle oldu?
disp(chess);

_chess = zeros(10,10);
[h,w] = size(_chess);
_chess(1:2:h, 2:2:w+1) = 1;
_chess(2:2:h+1, 1:2:w) = 1

%% Görüntü Okuma
% Üç Kanallı imge
I = imread("lena.png")

R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);


%GrayScale
Igray = rgb2gray(I);

figure;
imshow(I);
title('Original');

%Inew = (I*contrast) + parlaklık;

figure;
Inew = Igray * 2;
%imshow()

imshow(Igray - 100);
title('Siyah Beyaz');

%% Threshold Uygulama
figure;
imshow(Igray > 128,[]);
title('Threshold Uygulama');

figure;
%1 Transpoze alma
subplot(1,2,1);
I_transpose = Igray';
imshow(I_transpose);

%2 45 Derece Döndürme
subplot(1,2,2);
I45 = imrotate(I,45);
imshow(I45);

%% Mirror - 1
[h,w] = size(Igray);
new_gray = zeros(h,w);
for i=1:1:h
    for j=1:1:w
            new_gray(i, MOD(-j, w+1)) = Igray(i,j);
    end
end

figure;
imshow(uint8(new_gray));


%% Mirror - 2
new_gray_basic = zeros(h,w);
new_gray_basic = Igray(:,w:-1:1);
figure;
imshow(uint8(new_gray_basic));
title('Basic');

cYumusat = ones(5, 5) / 25; 
gaussian = fspecial ("gaussian", 11, 2);
sobel_y = fspecial('sobel');
sobel_x = fspecial('sobel')';

im_yumusak_kernel = imfilter(Igray, cYumusat, 'same');
im_yumusak_gauss = imfilter(Igray, gaussian, 'same');
im_sobel_y = imfilter(Igray, sobel_y, 'same');
im_sobel_x = imfilter(Igray, sobel_x, 'same');

figure;
subplot(2,2,1);
imshow(im_yumusak_kernel);
subplot(2,2,2);
imshow(im_yumusak_gauss);
subplot(2,2,3);
imshow(im_sobel_y);
subplot(2,2,4);
imshow(im_sobel_x);

pause;

%% 
  %%FFT
%%

%% FFT
fs = 100;               % sampling frequency
t = 0:(1/fs):(10-1/fs); % time vector
S = cos(2*pi*15*t);
n = length(S);
X = fft(S);
f = (0:n-1)*(fs/n);     %frequency range
power = abs(X).^2/n;    %power
pause(1);
figure;
plot(f,power)


%% FFT Shift
pause(1);
figure;
Y = fftshift(X);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
plot(fshift,powershift)

%% İmgelerde FFT Kullanımı
pause(1);
figure;
fft_image = fft2(Igray);
fft_shift_image = fftshift(fft_image);
subplot(2,1,1);
imshow(real(fft_image));
subplot(2,1,2);
imshow(real(fft_shift_image));

