% Testing Voigt Profile

clear
set_parameters
% cd minFunc_2012
% addpath(genpath(pwd))
% mexAll
% cd ..
% mex voigt0.c -lcerf
% mex voigt1.c -lcerf
dw = 0.25;
l1 =1548;
l2= 1550;
dl=200;
w1 = l1-dl:dw:l1+dl;
w2 = l2-dl:dw:l2+dl;
width=3;
padded_wavelengths1 = [linspace(min(w1), min(w1)+width*dw, width)';...
                      w1'; linspace(max(w1), l1+width*dw, width)']; 
padded_wavelengths2 = [linspace(min(w2), min(w2)+width*dw, width)';...
                      w2'; linspace(max(w2), l2+width*dw, width)'];
% z=0;
% for N=logspace(13,16,7)
%     clf()
%     ylim([-0.1 1])
%     title(sprintf('N:%.2e',N))
%     i=0;
%     c = ['b', 'r', 'g', 'k'];
%     for sigma=linspace(5e5,30e5,4)
%         i=i+1;
%         L1 = voigt0(padded_wavelengths1, z,N, 1, sigma);
%         L2 = voigt1(padded_wavelengths2, z, N, 1, sigma);
%         plot(w1, L1, 'Color',c(i))
%         hold on 
%         plot(w2, L2, 'Color',c(i))
%         hold on 
        
%         label{i} = sprintf('W1:%.4f, W2:%.4f, S:%.1e', trapz(w1, 1-L1), trapz(w2, 1-L2), sigma);
%     end
%     legend(label)
%     hold off
%     saveas(gca, sprintf('N-%.2e.png',N))
% end
L1 = voigt0(observed_wavelengths(padded_wavelengths1, 2.22),...
 1.85,10^14.87, 1, 50.80e5);
L2 = voigt1(observed_wavelengths(padded_wavelengths2, 2.22),...
 1.85, 10^14.87, 1, 50.80e5);
 figure()
%  plot(w1, L1)
%  hold on 
%  plot(w2, L2)
%  hold on 
 plot(w1, L1.*L2)
ylim([-0.1,1.1])