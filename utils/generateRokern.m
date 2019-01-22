 function rokern = generateRokern(N) 
        
        N= N*1.5;
        
        delta = [1.0, 0.0];
        k_not = [0.0, 0.0, 0.0];
        DCF   = [1.0];
        
        rokern = grid3_MAT(delta',k_not',DCF,N,1);
        
        % ROLLOFF
        % change to complex, shift, then fft
        rokern = squeeze(rokern(1,:,:,:) + 1j*rokern(2,:,:,:));
        rokern = fftn(rokern);
        rokern = fftshift(rokern,1);
        rokern = fftshift(rokern,2);
        rokern = fftshift(rokern,3);
        rokern = abs(rokern);
    end 