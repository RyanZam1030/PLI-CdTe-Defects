
% Clears command window
clc; 

%--------------------Initializing Photolumineces Data---------------------%
% Reading the text file. 
A = importdata('z400_7.txt','\t',1); 
% Grabs data from 2nd column of text file 
PLdata = A.data(:,2)'; 
% Smooths and flips the data 
PLdata = wrev(smooth(PLdata, 40, 'sgolay'));
% Converting the measurement to pixel. How many meters in a pixel
Xdata = (2.0000e-08 * (0:length(PLdata)-1)); 
% Finding the minimum value and column accociated with it
[M,I] = min(PLdata); 
Xdata(I);    
% Shifts the plot to match the min of the data.
Xdata = Xdata - Xdata(I); 

%-----------------Initializing Variables for Fucnction--------------------%
% Lifetime
tau_p = 1e-9;
% Resolution 
resolution = 0.25e-6  ;
% Total number of pixels in simulation
N = 100;   
% Number of pixels to scan/calculate 
Npix = N/2; 
% Reading % of laser power from text file
perc = str2double(A.colheaders(2)) / 100; 
% Power could be 20e-6 but is probably 60e-6 watts * 75%, 25%, 7%, or 2% 
power = 60e-6 * perc;
% Sample doping
Nd = 1e13;
% Leave 0
Nt = 0;
% Leave 0
Et = 0;

% Leave true
include_recycling_as_uniform_generation = true;
finalfignum = 15;
%-----------------------------Function Call-------------------------------% 
PLI = cPL_Relax(tau_p, resolution, N, power, Nd, Nt, Et , Npix, true, finalfignum);
% Modeled value for Xdata and reverses vector
X = wrev((0:Npix - 2) * resolution);  
figure(finalfignum);
plot(X, PLI / max(PLI), Xdata, PLdata/max(PLdata));

%-----------------Function Photoluminescent Intensity(PLI)----------------%
function PLI = cPL_Relax(tau_p, resolution, pixlesInSim, power,... 
                         sampleDope, Nt, Et , PixelsToScan,... 
                         include_recycling_as_uniform_generation, finalFigNum)
fram = 0;                                                                  %Initializing 
F = getframe;                                                              %Initializing 
filmThick = 2e-6;                                                          %film thickness
Nt = Nt / filmThick;                                                       %value is 0 
tau_no = 1 / (Nt * 4e-20 * 3.8e5);                                         %value is inf
tau_po = 1 / (Nt * 4e-19 * 1.3e5);                                         %value is inf
maxItCount = 1e6;                                                          %maximum number of allowed iteration
sampleDope = sampleDope * maxItCount; 
Nc = 7.3e17 * maxItCount;
Nv = 1.8e19 * maxItCount;
Eg = 1.5 / 0.026;
Et = Et / 0.025;                                                           %value is 0
n1 = Nc * exp(-(Eg - Et));  
p1 = Nv * exp(-(Et));
brad = 1e-16;
include_recycling_as_diffusion = false; 
twinoffset = 0;                                                            %set to 0 for simple dislocation
U = zeros(pixlesInSim, pixlesInSim);                                       %Initializing a matrix
dist = @(di, dj) sqrt((di).^2 + (dj).^2) * resolution;                     %distance from center
for i = 1:pixlesInSim
  
   for j = 1:pixlesInSim
      U(i,j) = dist(i - pixlesInSim / 2, j - pixlesInSim / 2);
   end
   
end

spot_rad = 1 * 0.25e-6;                                                    %spot radius
Gnorm = exp(-U.^2 / (2 * spot_rad^2)) / (2 * pi * spot_rad^2);             %2-D Gaussian distribution (m-2). Integrates to 1 over area (Ndx)^2;
Generatoin = power * Gnorm;                                                %Generation: W/m2
Generatoin = Generatoin / ((1240 / 488) * 1.602e-19);                      %carrier/m2 per s at 488 nm (or 514 nm)
Generatoin = Generatoin / (filmThick);                                     %Uniform generation with depth: carrier/m3 per s.

%------------------Initializing boundary conditions-----------------------%
U = zeros(pixlesInSim, pixlesInSim);                                       %Initializing a matrix
Generatoin(1,:) = 0; Generatoin(pixlesInSim,:) = 0;  
Generatoin(:,1) = 0; Generatoin(:,pixlesInSim) = 0;
Unew = U;                                                                  %Carrier concentration 
Xpos = round((pixlesInSim - 1) / 2);                                       %Xpos, Ypos are coordinates of dislocation. Laser spot is center.
PLI = ones(1, Xpos - 1);                                                   %Photo luminecent intensity (list for results)
X = (0:Xpos - 2) * resolution;  
X = wrev(X);                                                               %reverses vector

for Ypos = 2:PixelsToScan                                                  %Loop over dislocation position
   
   Twinpos = (Ypos - round(twinoffset / resolution));
   if(Twinpos < 1)
      Twinpos = Ypos;
   end
   U = circshift(U, 1, 1);
   U(1,:) = 0; U(pixlesInSim,:) = 0; U(:,1) = 0; U(:,pixlesInSim) = 0;
   U(Ypos, Xpos) = 0;
   U(Twinpos, Xpos) = 0;
   Ncount = 0;
   loop = true;
   while loop                                                              %Loop until tolerance reached
      
      Ncount = Ncount + 1;
      mobility = (1e-4 * 800 * 60 *(U + U + sampleDope))./... 
                 (800 * (U + sampleDope) + 60 * U);
      D = 0.026 * mobility;
      if(Nt > 0)
         B = brad + 1./ ((U + p1) * tau_no + (U + sampleDope + n1) * tau_po);
      else
         B = brad;
      end
      
      if(include_recycling_as_diffusion)
         recycleAbs = 1e5;
         Dtot = D + (pi / (4 * recycleAbs))^2 * brad * (U + sampleDope)/ 2;
         B = B - brad * (1 - 0.2);
      else
         Dtot = D;
      end
      
      if(include_recycling_as_uniform_generation)
         Rrad = U.* (U + sampleDope) * brad;
         Gtot = Generatoin + (1 - 0.2) * sum(sum(Rrad)) / pixlesInSim^2;
      else
         Gtot = Generatoin;
      end
      
      gammB = Dtot./ (resolution^2 * B);
      Unear = (U(Ypos + 1, Xpos) + U(Ypos - 1, Xpos) + ...
               U(Ypos, Xpos + 1) + U(Ypos, Xpos - 1));
      rate = exp(-1e5 * sqrt(8.85e-12 * 0.025 / ...
                (1.602e-19 * (Unear / 4 + sampleDope)))) / (tau_p);        %sub-pixel scale screening
      Unew = -(+sampleDope + 4 * gammB) / 2 + ...
               sqrt((sampleDope + 4 * gammB).^2 / 4 + Gtot./B + ...
               gammB.* (circshift(U, 1, 1) + circshift(U, -1, 1) + ...
               circshift(U, 1, 2) + circshift(U, -1, 2)));
      Unew(1,:) = 0; 
      Unew(pixlesInSim,:) = 0;  
      Unew(:,1) = 0;  
      Unew(:,pixlesInSim) = 0;
      if(tau_p == 0)
         Ndark = 0;
         Unew(Ypos - Ndark:Ypos + Ndark, Xpos - Ndark:Xpos + Ndark) = 0;
         Unew(Ypos,Xpos) = 0;
         Unew(Twinpos, Xpos) = 0;
      else
         gammBi = gammB(Ypos, Xpos);
         Utau = -(rate + 4*B*gammBi - ...
             sqrt(B^2*sampleDope^2 + 8*B^2*sampleDope*gammBi + 16*B^2*gammBi^2 ...
             + 4*Unear*B^2*gammBi + 2*B*sampleDope*rate + 8*B*gammBi*rate + ...
             4*Gtot(Ypos, Xpos)*B + rate^2) + B*sampleDope)/(2*B);
         Unew(Ypos, Xpos) = Utau;
      end
      
      ERR =  sum(sum(abs(Unew - U))) / sum(sum((Unew)));
      U = Unew;
      if(ERR < 1e-6)                                                       %rate of change in terms of fractional change
         loop = false;
         disp(['solution reaches steady state in ',num2str(Ncount) ,' steps'])
      end
      
      if (mod(Ncount, maxItCount / 20)==0 || ~loop)                        %displays movie frame every 100 time steps
         
         if(Ncount >= maxItCount)
            disp(['solution does not reach steady state in ',num2str(maxItCount),' steps'])
         end
         
         dgh = graphit(U);                                                 %Function Call
      end
      
   end
   pause on;
   fram = fram + 1;
   F(fram) = getframe(dgh);
   PLI(Ypos-1) = brad * filmThick * resolution.^2.* sum(sum(Gnorm.* U.^2));
   figure(finalFigNum);
   plot(X, PLI / max(PLI));
   pause (0.01);
   pause off;
end

%--------------------Display a movie of heat diffusion--------------------%
v = VideoWriter('RadRecombDiff.avi');
open(v);
writeVideo(v, F);
close(v);
end

%----------------------Function diffusion graph(dgh)-----------------------%
function dgh = graphit(U)
pause on;
N = length(U);
dgh = figure(1);
Ugraph = horzcat(flip(U,2),U);
surf(Ugraph);
axis([1 N 1 N ]);
h = gca;
get(h,'FontSize') ;
set(h,'FontSize',12)
colorbar('location','eastoutside','fontsize',12);
xlabel('X','fontSize',12);
ylabel('Y','fontSize',12);
title('Diffusion','fontsize',12);
set(dgh, 'color', 'white');
pause(0.01);
pause off;
end

%----------------------------------Notes----------------------------------%
% IP over 1e-4 Watts is radiative-recombination-dominated when 
% limited to a 15 micron radius region. z404 can be replaced, 
% but z403 only has one huge spot image.


