%% EXPERIMENT 23
% DSP-Based Image Enhancement, Compression and Edge Detection Using MATLAB
% Input: cameraman.tif (256 x 256 grayscale sample image)
% Processing: median filtering, contrast stretching, JPEG compression, Canny edges
% Required toolbox: Image Processing Toolbox

clc;
clear;
close all;

%% 1. INPUT PARAMETERS
imgName = 'cameraman.tif';
medSize = [3 3];
jpegQuality = 50;
cannyTh = [0.1 0.2];
MAXI = 255;

%% 2. LOAD INPUT IMAGE (ROBUST SEARCH)
% Search the current folder and MATLAB path first.
imgPath = '';

if isfile(imgName)
    imgPath = imgName;
else
    foundPath = which(imgName);
    if ~isempty(foundPath) && isfile(foundPath)
        imgPath = foundPath;
    end
end

% Try MATLAB's standard sample-image location, if available.
if isempty(imgPath)
    candidate = fullfile(matlabroot, 'toolbox', 'images', 'imdata', imgName);
    if isfile(candidate)
        imgPath = candidate;
    end
end

% If not found, stop with a useful message before later calculations.
if isempty(imgPath)
    error(['Cannot find %s. Upload cameraman.tif to the same folder as this script ' ...
           'in MATLAB Drive, then run the script again.'], imgName);
end

I = imread(imgPath);

if size(I, 3) == 3
    I = rgb2gray(I);
end
if ~isa(I, 'uint8')
    I = im2uint8(I);
end

[rows, cols] = size(I);
fprintf('============================================================\n');
fprintf(' EXPERIMENT 23: IMAGE ENHANCEMENT, COMPRESSION AND EDGE DETECTION\n');
fprintf('============================================================\n');
fprintf('Image file        : %s\n', imgPath);
fprintf('Image size        : %d x %d\n', rows, cols);
fprintf('Image class       : %s\n', class(I));
fprintf('Minimum intensity : %d\n', min(I(:)));
fprintf('Maximum intensity : %d\n', max(I(:)));
fprintf('Mean intensity    : %.2f\n', mean(double(I(:))));
fprintf('Std. deviation    : %.2f\n', std(double(I(:))));

%% 3. MEDIAN FILTERING AND CONTRAST ENHANCEMENT
I_med = medfilt2(I, medSize);
I_enh = imadjust(I_med, stretchlim(I_med), []);

%% 4. JPEG COMPRESSION AND RECONSTRUCTION
jpgFile = fullfile(tempdir, 'cameraman_q50.jpg');
imwrite(I, jpgFile, 'jpg', 'Quality', jpegQuality);
I_rec = imread(jpgFile);

origBytes = numel(I);  % 8-bit grayscale image: 1 byte per pixel
jpgInfo = dir(jpgFile);
compBytes = jpgInfo.bytes;
CR = origBytes / compBytes;
bpp = 8 * compBytes / numel(I);

%% 5. CANNY EDGE DETECTION
E_orig = edge(I, 'canny', cannyTh);
E_enh  = edge(I_enh, 'canny', cannyTh);
E_rec  = edge(I_rec, 'canny', cannyTh);

%% 6. IMAGE QUALITY METRICS
% PSNR is infinite when the images are identical (MSE = 0).
MSE = @(A,B) mean((double(A(:)) - double(B(:))).^2);
PSNR = @(A,B) localPSNR(MSE(A,B), MAXI);

MSE_med = MSE(I, I_med);
PSNR_med = PSNR(I, I_med);

MSE_enh = MSE(I, I_enh);
PSNR_enh = PSNR(I, I_enh);

MSE_rec = MSE(I, I_rec);
PSNR_rec = PSNR(I, I_rec);

NIE_med = MSE_med / MAXI^2;
NIE_enh = MSE_enh / MAXI^2;
NIE_rec = MSE_rec / MAXI^2;

SSIM_rec = ssim(I_rec, I);
edgePix = [nnz(E_orig), nnz(E_enh), nnz(E_rec)];

%% 7. FIGURE 23A: IMAGE PROCESSING RESULTS
figure('Name', 'Experiment 23A: Image Processing Results', ...
    'Color', 'w', 'Position', [40 40 1200 750]);

subplot(2,3,1);
imshow(I);
title('Original image');

subplot(2,3,2);
imshow(I_med);
title('After 3 x 3 median filter');

subplot(2,3,3);
imshow(I_enh);
title('Enhanced (median + imadjust)');

subplot(2,3,4);
imshow(I_rec);
title(sprintf('JPEG Q=%d reconstructed (CR = %.2f)', jpegQuality, CR));

subplot(2,3,5);
imshow(E_orig);
title('Canny edge map (original)');

subplot(2,3,6);
imshow(E_rec);
title('Canny edge map (JPEG reconstructed)');

%% 8. FIGURE 23B: HISTOGRAMS
figure('Name', 'Experiment 23B: Histograms', ...
    'Color', 'w', 'Position', [80 80 900 350]);

subplot(1,2,1);
imhist(I);
title('Histogram - original');

subplot(1,2,2);
imhist(I_enh);
title('Histogram - enhanced');

%% 9. FIGURE 23C: COMPRESSION ERROR AND ENHANCED EDGES
figure('Name', 'Experiment 23C: Compression Error and Enhanced Edges', ...
    'Color', 'w', 'Position', [120 120 900 380]);

subplot(1,2,1);
imshow(abs(double(I) - double(I_rec)), []);
colorbar;
title('|Original - JPEG reconstructed|');

subplot(1,2,2);
imshow(E_enh);
title('Canny edge map (enhanced image)');

%% 10. FIGURE 23D: JPEG RATE-DISTORTION ANALYSIS
qList = 10:10:90;
psnrQ = zeros(size(qList));
crQ = zeros(size(qList));

for k = 1:numel(qList)
    tempJpg = fullfile(tempdir, sprintf('cameraman_q%d.jpg', qList(k)));
    imwrite(I, tempJpg, 'jpg', 'Quality', qList(k));

    I_temp = imread(tempJpg);
    tempInfo = dir(tempJpg);

    psnrQ(k) = PSNR(I, I_temp);
    crQ(k) = origBytes / tempInfo.bytes;
end

figure('Name', 'Experiment 23D: JPEG Rate-Distortion', ...
    'Color', 'w', 'Position', [160 160 900 380]);

subplot(1,2,1);
plot(qList, psnrQ, '-o', 'LineWidth', 1.4);
grid on;
xlabel('JPEG quality');
ylabel('PSNR (dB)');
title('PSNR vs JPEG quality');

subplot(1,2,2);
plot(qList, crQ, '-s', 'LineWidth', 1.4);
grid on;
xlabel('JPEG quality');
ylabel('Compression ratio');
title('Compression ratio vs JPEG quality');

%% 11. OUTPUT PARAMETER TABLE
fprintf('\n================ OUTPUT PARAMETER TABLE ================\n');
fprintf('%-28s %12s %12s %14s\n', ...
    'Image', 'MSE', 'PSNR (dB)', 'MSE/MAXI^2');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%-28s %12.2f %12.2f %14.6f\n', ...
    'Median filtered', MSE_med, PSNR_med, NIE_med);
fprintf('%-28s %12.2f %12.2f %14.6f\n', ...
    'Enhanced (imadjust)', MSE_enh, PSNR_enh, NIE_enh);
fprintf('%-28s %12.2f %12.2f %14.6f\n', ...
    sprintf('JPEG Q=%d reconstructed', jpegQuality), MSE_rec, PSNR_rec, NIE_rec);

fprintf('\n--- Compression results (JPEG Q = %d) ---\n', jpegQuality);
fprintf('Original size      : %d bytes\n', origBytes);
fprintf('Compressed size    : %d bytes\n', compBytes);
fprintf('Compression ratio  : %.2f : 1\n', CR);
fprintf('Bits per pixel     : %.3f bpp\n', bpp);
fprintf('SSIM (JPEG vs original): %.4f\n', SSIM_rec);

fprintf('\n--- Edge detection (Canny [%.1f %.1f]) ---\n', cannyTh(1), cannyTh(2));
fprintf('Edge pixels - original image        : %d\n', edgePix(1));
fprintf('Edge pixels - enhanced image        : %d\n', edgePix(2));
fprintf('Edge pixels - JPEG reconstructed    : %d\n', edgePix(3));

fprintf('\n--- JPEG rate-distortion results ---\n');
fprintf('%12s %14s %20s\n', 'Quality', 'PSNR (dB)', 'Compression ratio');
for k = 1:numel(qList)
    fprintf('%12d %14.2f %20.2f\n', qList(k), psnrQ(k), crQ(k));
end

fprintf('\nRESULT: Image enhancement, JPEG compression and edge detection completed.\n');
fprintf('============================================================\n');

%% LOCAL FUNCTION
function value = localPSNR(mseValue, peakValue)
    if mseValue == 0
        value = Inf;
    else
        value = 10 * log10((peakValue^2) / mseValue);
    end
end
