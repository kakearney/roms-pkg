% Demo source code for interp1dcol
clc;
clear all;
close all;
X = 0 : 2;
nCols = 10000;
Xarr = repmat(X', [1 nCols]);

Varr = zeros(length(X), nCols);
for colIdx =  1: nCols
    Varr(:, colIdx) = sin(X' +  colIdx / nCols * 2 * pi);
end

Xq = 0 : 0.25 : 2;
Xqarr = repmat(Xq', [1 nCols]);
Vq = zeros(length(Xq), nCols);
tic;
for colIdx = 1 : nCols
    Vq(:, colIdx) = interp1q(Xarr(:, colIdx), Varr(:, colIdx), Xqarr(:, colIdx));
end
toc;

% Now, call the interp1dCol
tic;
Vq2 = fastinterpcol(Xarr, Varr, Xqarr);
toc;
disp(['Error between interp1q and interpcol: ' num2str(norm(Vq - Vq2, 'fro'))]);
figure(1);
plot(Xq, Vq(:, 1));
hold on;
plot(Xq, Vq2(:, 1), 'r');
legend('Interp1q', 'Interpcol');
