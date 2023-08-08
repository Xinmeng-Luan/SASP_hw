%------------------------------------------%
%        *** SSSP - HOMEWORK #3 ***        %
%------------------------------------------%
%     Emulation of the Leslie Speaker      %
%------------------------------------------%
% Name: Iaccarino - Luan                   %
% Student ID: 10868500 - 10876787          %
%------------------------------------------%

function [y,y_lpf,y_hpf,y_hp_sdf] = leslie_o(x, Fs, freq)
%Leslie Speaker Emulation
%
% J. Pekonen et al. Computationally Efficient Hammond Organ Synthesis
% in Proc. of the 14th International Conference on Digital Audio
% Effects(DAFx-11), Paris, France, Sept. 19-23, 2011

% length of the input signal
N = length(x);
N_channel = size(x,2);
% time axis
t = 0: 1/Fs: (N/Fs);
t = t';
% 
% pad_t = zeros(N_sdf_t, size(x,2));
% x_padded_t = [pad_t x];
% pad_b = zeros(N_sdf_b, size(x,2));
% x_padded_b = [pad_b x];

% global modulator parameters
alpha=0.9;
% tremble spectral delay filter parameter 
Ms_t=0.2;
Mb_t=-0.75;
N_sdf_t=4;
% bass spectral delay filter parameter 
Ms_b=0.04;
Mb_b=-0.92;
N_sdf_b=3;

% cross-over network design
fc=800;                 % cutoff frequency
N_cn = 4;
%TODO: compute the coefficients for the two 4th order butterworth filters
%with cutoff frequency fc
[b_lp, a_lp]= butter(N_cn, fc/(Fs/2), "low"); %LPF design
[b_hp, a_hp]= butter(N_cn, fc/(Fs/2), "high");  %HPF design

% allocate input and output buffers for IIR filters
% hp filter buffers
hpf.state=zeros(N_cn,N_channel);    %output (final state of the filter)
hpf.in=zeros(N_cn+1,N_channel);       %input (initial conditions of the filter)
% lp filter buffers
lpf.state=zeros(N_cn,N_channel);
lpf.in=zeros(N_cn+1,N_channel);
% treble sdf filter buffers
sdf_h.state=zeros(N_sdf_t+1,N_channel);
sdf_h.in=zeros(N_sdf_t+1,N_channel);
% bass sdf filter buffers
sdf_b.state=zeros(N_sdf_b+1,N_channel);
sdf_b.in=zeros(N_sdf_b+1,N_channel);

% modulators
m_sin_b = sin(2*pi*freq.*t);
m_sin_t = sin(2*pi*(freq+0.1).*t);
m_b = Ms_b .* m_sin_b + Mb_b; % bass modulator
m_t = Ms_t .* m_sin_t + Mb_t; % tremble modulator


y_lpf = zeros(N, N_channel);
y_lpf_sdf = zeros(N, N_channel);
y_hpf = zeros(N, N_channel);
y_hpf_sdf = zeros(N, N_channel);
y_lp_am = zeros(N, N_channel);
y_hp_am = zeros(N, N_channel);

y = zeros(N, N_channel);

%sample processing
for n=1:N

    % compute crossover network filters outputs
    %% y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) + b3*x(n-3) + b4*x(n-4) -
    %%      - a1*y(n-1) - a2*y(n-2) - a3*y(n-3) - a4*y(n-4)   
    
    %% lowpass butterworth %%
    %  x(n-delay)
    for j = 1:N_cn+1
        if n-(j-1)>0
            lpf.in(j,:) = x(n - (j-1),:);
        end
    end
    
    % input related terms
    butterSumX = zeros(1, N_channel);
    for d = 1:N_cn+1
        butterSumX = butterSumX + b_lp(d).*lpf.in(d,:);
    end

    % input related terms
    butterSumY = zeros(1, N_channel);
    for d = 1:N_cn
        butterSumY = butterSumY - a_lp(d+1).*lpf.state(d,:);
    end

    %computation of y_lpf
    y_lpf(n,:) = butterSumX + butterSumY;   

    
    %  y(n-delay)
    % index j corresponds to y(n-(j-1)) that becomes y(n-j) 
    % in the following iteration of the processing loop
    for j = 1:N_cn
        if n-(j-1)>0
            lpf.state(j,:) = y_lpf(n - (j-1),:);
        end
    end
    

    %% highpass butterworth %%
    %  x(n-delay)
    for j = 1:N_cn+1
        if n-(j-1)>0
            hpf.in(j,:) = x(n - (j-1),:);
        end
    end
    
    % input related terms
    butterSumX = zeros(1, N_channel);
    for d = 1:N_cn+1
        butterSumX = butterSumX + b_hp(d).*hpf.in(d,:);
    end

    % input related terms
    butterSumY = zeros(1, N_channel);
    for d = 1:N_cn
        butterSumY = butterSumY - a_hp(d+1).*hpf.state(d,:);
    end

    %computation of y_lpf
    y_hpf(n,:) = butterSumX + butterSumY;   

    
    %  y(n-delay)
    % index j corresponds to y(n-(j-1)) that becomes y(n-j) 
    % in the following iteration of the processing loop
    for j = 1:N_cn
        if n-(j-1)>0
            hpf.state(j,:) = y_hpf(n - (j-1),:);
        end
    end
    
%     [y_lpf, lpf.state] = filter(b_lp, a_lp, x(n));
    


    %% compute bass SDF output %%

    % putting the last N_sdf_t samples into the buffer
    % the sdf.in buffer corresponds to the x signal of the difference eq.
%     for j = 1:N_sdf_b+1
%         if n-(j-1)>0
%             sdf_b.in(j,:) = y_lpf(n-(j-1),:);
%         end
%     end                         
% 

    sdf_b.in(2:end,:) = sdf_b.in(1:end-1,:);
    sdf_b.in(1,:) = y_lpf(n);
                                          
    % summation from 1 to N+1 -> corresponds to the summation from 0 to N
    % of the analytical formulation of the SDFs
    temp_sum = zeros(1,N_channel);
    for i = 1:N_sdf_b+1
        bin_coeff = nchoosek(N_sdf_b,i-1);       %i-1=i of the pdf
        temp_sum = bin_coeff * m_b(n+1)^(i-1) * ( sdf_b.in(N_sdf_b+2-i,:) - sdf_b.state(i,:) );
        %the final value gets updated each cycle
        y_lpf_sdf(n,:) = y_lpf_sdf(n,:) + temp_sum;
    end

    % copying the last values of y into the buffer
%     for j = 1:N_sdf_b+1
%         if n-(j-1)>1 && (j-1)>0
%             sdf_b.state(j,:) = y_lpf_sdf(n-(j-1),:);
%         end
%     end

    sdf_b.state(3:end,:) = sdf_b.state(2:end-1,:);
    sdf_b.state(2,:) = y_lpf_sdf(n);

    %% compute treble SDF output
    
    % putting the last N_sdf_t samples into the buffer
    % the sdf.in buffer corresponds to the x signal of the difference eq.
                       
    sdf_h.in(2:end,:) = sdf_h.in(1:end-1,:);
    sdf_h.in(1,:) = y_hpf(n);                              

    % summation from 1 to N+1 -> corresponds to the summation from 0 to N
    % of the analytical formulation of the SDFs
    temp_sum = zeros(1,N_channel);
    for i = 1:N_sdf_t+1
        bin_coeff = nchoosek(N_sdf_t,i-1);       %i-1=i of the pdf
        temp_sum = bin_coeff * m_t(n+1)^(i-1) * ( sdf_h.in(N_sdf_t+2-i,:) - sdf_h.state(i,:) );
        %the final value gets updated each cycle
        y_hpf_sdf(n,:) = y_hpf_sdf(n,:) + temp_sum;
    end

    %copying the last values of y into the buffer
    sdf_h.state(3:end,:) = sdf_h.state(2:end-1,:);
    sdf_h.state(2,:) = y_hpf_sdf(n);

    %% implement AM modulation blocks
    y_lp_am(n,:) = ( 1 + alpha*m_b(n+1) ) * y_lpf_sdf(n,:);
    y_hp_am(n,:) = ( 1 + alpha*m_t(n+1) ) * y_hpf_sdf(n,:);

    y(n,:) = y_lp_am(n,:) + y_hp_am(n,:);

end

end


