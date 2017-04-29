close all; clear; clc; warning('off','all'); tic

% Inisiasi awal
d = DST_class;                      % Referensikan file class OOP
d.equation = 'parabolic';           % Pilih equation, 'parabolic' atau 'elliptic'
d.matriks = csvread('DST_N1.csv');  % Baca matriks pembobot dari file eksternal

% Input semua parameter
d.N = 5;                % Perhalusan
d.dt = 0.1;             % Delta time
d.rho = 2702;           % Massa jenis
d.cp = 0.903;           %
d.k = 20;               %
d.fluks = 7900;         % Fluks kalor
d.TLing = 25;           % Suhu lingkungan
d.panjang = 0.5;        % Dimensi plat
d.lebar = 0.2;
d.tebal = 0.1;
d.temp_A = [3,3];       % Titik A
d.temp_B = [2,10];      % Titik B
d.fluks_kiri = [3,7];   % Titik C
d.fluks_kanan = [3,11]; % Titik D

% Olah data
d.matriks = d.matriks';
d = find_SqrSmooth(d);      % Memperbesar matriks pembobot
[d.n,d.m] = size(d.Sqr);    % Ukuran matriks pembobot
d.dx = d.panjang/(d.n-1);   % Panjang per segmen
d = find_Q(d);              % Menentukan fluks kalor tiap segmen
d = find_Fo(d);             % Menentukan Fo untuk transien

% Menentukan jumlah iterasi transien
% Mendefinisikan suhu awal transien
iter = 1;
T{1} = d.TLing*ones(d.m,d.n);
if strcmp(d.equation,'parabolic')
    iter = 100;
end

% Mencari matriks A
A = cell(d.m);
[A{:}] = deal(zeros(d.n));
for i = 1:d.m
    if i > 1
        A{i,i-1} = find_side(d,i);
    end
    A{i,i} = find_A(d,i);
    if i < d.m
        A{i,i+1} = find_side(d,i);
    end
end
A = cell2mat(A);

% Looping transien
for k = 1:iter
    % Mencari vektor C
    C = cell(d.m,1);
    [C{:}] = deal(zeros(d.n,1));
    for i = 1:d.m
        C{i} = find_C(d,i,T{k});
    end
    C = cell2mat(C);
    
    % Hasil simulasi
    T{k+1} = find_T(d,A,C);             % Mencari matriks T
    plot_Mesh(d,T{k+1},k);              % Memplot grafiknya
    d = find_temperaturAB(d,T{k+1},k);  % Temperatur di A dan B
    d = find_fluksCD(d,T{k+1},k);       % Fluks di CD
end

clear fluks_kiri fluks_kanan i k iter matriks; toc