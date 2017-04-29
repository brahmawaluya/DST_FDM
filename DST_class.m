classdef DST_class
    properties
        equation
        N
        dt
        rho
        cp
        alpha
        Fo
        k
        fluks
        TLing
        panjang
        lebar
        tebal
        temp_A
        temp_B
        fluks_kiri
        fluks_kanan
        matriks
        n
        m
        dx
        Q
        Sqr
        T_A
        T_B
        fluksCD
    end
    methods
        
        % Menentukan fluks kalor tiap segmen
        function d = find_Q(d)
            Q_t = d.fluks/(d.tebal*d.dx);
            d.Q = 2*Q_t*d.dx/(d.k*d.N);
        end
        
        % Menentukan Fo untuk transien
        function d = find_Fo(d)
            d.alpha = d.k/(d.rho*d.cp);
            d.Fo = d.alpha*d.dt/d.dx^2;
        end
        
        % Memperbesar matriks pembobot
        function d = find_SqrSmooth(d)
            % Membuat matriks baru yang telah diperbesar (d.Sqr)
            [p,l] = size(d.matriks);
            d.Sqr = zeros(p*d.N-d.N+1,l*d.N-d.N+1);
            % Petakan semua nilai ke lokasi (x y) di matriks yang baru
            for j = 1:l
                for i = 1:p
                    x = d.N*(i-1)+1;
                    y = d.N*(j-1)+1;
                    % Assign nilai-nilai di matriks baru tsb
                    switch d.matriks(i,j)
                        % Nilai 1 akan menyebar ke semua arah
                        case 1,
                            for a = -(d.N-1):d.N-1
                                for b = -(d.N-1):d.N-1
                                    d.Sqr(x+a,y+b) = d.matriks(i,j);
                                end
                            end
                        % Nilai 2 dan 3 menyebar ke kiri kanan/atas bawah    
                        case {2,3},
                            if j == 1 || j == l;
                                for a = -(d.N-1):d.N-1
                                    d.Sqr(x+a,y) = d.matriks(i,j);
                                end
                            elseif i == 1 || i == p;
                                for b = -(d.N-1):d.N-1
                                    d.Sqr(x,y+b) = d.matriks(i,j);
                                end
                            end
                        % Nilai 6 tidak menyebar    
                        case 6,
                            d.Sqr(x,y) = d.matriks(i,j);
                    end
                    
                    % Mencari titik A dan B pada matriks yang diperbesar
                    if  j == d.temp_A(1) && i == d.temp_A(2)
                        d.temp_A(1) = y;
                        d.temp_A(2) = x;
                    elseif  j == d.temp_B(1) && i == d.temp_B(2)
                        d.temp_B(1) = y;
                        d.temp_B(2) = x;
                    end
                    
                    % Mencari titik C dan D pada matriks yang diperbesar
                    if j == d.fluks_kiri(1) && i == d.fluks_kiri(2)
                        d.fluks_kiri(1) = y;
                        d.fluks_kiri(2) = x;
                    elseif j == d.fluks_kanan(1) && i == d.fluks_kanan(2)
                        d.fluks_kanan(1) = y;
                        d.fluks_kanan(2) = x;
                    end   
                end
            end
        end
        
        % Mencari matriks A (diagonal utama)
        function A = find_A(d,i)
            A = zeros(d.n);
            for j = 1:d.n
                % Nilai default untuk lingkungan
                A(j,j) = -1;
                % Nilai untuk pembobot 1,2,3,6
                if d.Sqr(j,i) ~= 0
                    flux = 1;
                    if mod(d.Sqr(j,i),3) == 0
                        flux = 2;
                    end
                    
                    if j ~= 1
                        A(j,j-1) = flux;
                    end
                    if strcmp(d.equation,'parabolic')
                        A(j,j) = -(4 + 1/d.Fo);
                    else
                        A(j,j) = -4;
                    end
                    if j ~= d.n
                        A(j,j+1) = flux;
                    end
                end
            end      
        end
        
        % Mencari matriks A (diagonal sisi)
        function side = find_side(d,i)
            side = zeros(d.n);
            for j = 1:d.n
                if d.Sqr(j,i) ~= 0
                    ins = 1;
                    if mod(d.Sqr(j,i),2) == 0
                        ins = 2;
                    end
                    side(j,j) = ins;               
                end
            end
        end
        
        % Mencari vektor C
        function C = find_C(d,i,T)
            C = zeros(d.n,1);
            for j = 1:d.n
                switch d.Sqr(j,i)
                    case 0,
                        C(j) = -d.TLing;
                    case {1,2},
                        if strcmp(d.equation,'parabolic')
                            C(j) = -T(i,j)/d.Fo;
                        end
                    case {3,6},
                        if strcmp(d.equation,'parabolic')
                            C(j) = (-T(i,j) - (d.fluks * d.alpha * d.dt)/(d.k * d.dx))/d.Fo;
                        else
                            C(j) = -d.Q;
                        end
                end
            end
        end
        
        % Mencari matriks T
        function T = find_T(d,A,C)
            T = inv(A)*C;
            T = vec2mat(T,d.n);
        end
        
        % Memplot grafik simulasi
        function plot_Mesh(d,T,k)
            % Konfigurasi agar seluruh matriks T dapat dilihat di mesh
            T(d.m+1,d.n+1) = 0;
            T = circshift(T,1,1);
            
            % Buat mesh
            [x,y] = meshgrid(0:d.panjang/d.n:d.panjang, 0:d.lebar/d.m:d.lebar);
            
            % Konfigurasi untuk bagian plat yang kosong
            a = 4;
            b = 10;
            for i = (6/10*d.N*b)+3:(d.N*b)+2
                for j = 1:(2/4*d.N*a)
                    x(j,i) = NaN;
                    y(j,i) = NaN;
                end
            end
            
            % Plot setiap kali grafik selesai
            h = surface(x,y,flipud(T));
            set(h, 'edgecolor','none')
            title('Plot Temperatur Steady State 2D');
            colormap('hot'); colorbar;
            axis([0 d.panjang 0 d.lebar]);
            drawnow;
            %saveas(gcf,strcat('img_',num2str(k),'.jpg'));
        end
        
        % Temperatur di A dan B
        function d = find_temperaturAB(d,T,k)
            d.T_A(k) = T(d.temp_A(1),d.temp_A(2));
            d.T_B(k) = T(d.temp_B(1),d.temp_B(2));
            %fprintf('Temperatur A = %4.2f\nTemperatur B = %4.2f\n',T_A,T_B);
        end
        
        % Fluks di CD
        function d = find_fluksCD(d,T,k)
            fluxCD = 0;
            for i = d.fluks_kiri(2):d.fluks_kanan(2)
                fluxCD = fluxCD + d.k/d.dx*(T(d.fluks_kiri(1)-1,i)-T(d.fluks_kiri(1),i));
            end
            d.fluksCD(k) = fluxCD/d.N;
            %fprintf('Fluks di CD = %e\n',Q);
        end
        
    end
end