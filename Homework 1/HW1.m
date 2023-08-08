clc;
clear all;
close all;

%% Legend
% p: number past samples
% n: current time
% k: k_th coefficient
% M: segment length
% a: filter coefficients

% speech: 5ms, 

<<<<<<< Updated upstream
%%
% [hPiano, ePiano] = LPCFilter("piano.wav");
% [hSpeech, eSpeech] = LPCFilter("speech.wav");
[predict_piano] = LPCFilter("piano.wav");
predict_piano_reshape = reshape(predict_piano', 1, size(predict_piano,1)*size(predict_piano,2));
sound(predict_piano_reshape)
=======
% Import the files
[signal, fs] = audioread("piano.wav");

% 5 ms is taken from lesson as example segment length
M = floor(5e-3*fs); % How many samples in each segment


%%%%//// Method by Xinmeng
num_segment = ceil(length(signal)/M);
num_pad = num_segment* M -length(signal);
paddedSignal = padarray(signal,[num_pad 0],0,'post');
s = reshape(paddedSignal,M,num_segment)';
% 
pred =zeros(size(s));
err =zeros(size(s));

for ss = 1:num_segment
    [pred(ss,:) err(ss,:)] = lpc(s(ss,:), 219);
end

audioOut =zeros(size(s));
for ss = 1:num_segment
    audioOut(ss,:) = filter([0 -pred(ss,2:end)], 1, s(ss,:));
end

audioOut_reshape = reshape(audioOut',[848980 1]);

% sound(audioOut_reshape, 44100); 
>>>>>>> Stashed changes
%%
% [hSpeech, eSpeech] = LPCFilter("speech.wav");
% [a] = LPCFilter("piano.wav");

load a.mat;

a_exp1 = ones(size(a,1),1);
a_exp =[a_exp1 -1.*a];



H = zeros(size(a_exp));
for ii = 1:size(a_exp,1)
    H(ii,:) = freqz(1, a_exp(ii,:), 220);
end
A = 1./H;

a_exp0 = zeros(size(a,1),1);
a_exp_p =[a_exp0 a];

P = zeros(size(a_exp_p));
for ii = 1:size(a_exp_p,1)
    P(ii,:) = freqz(a_exp_p(ii,:),1,220);
end


%%%%%%%%

% Import the files
[signal, fs] = audioread("piano.wav");

% 5 ms is taken from lesson as example segment length
M = floor(5e-3*fs); % How many samples in each segment

%%%%//// Method by Xinmeng
num_segment = ceil(length(signal)/M);
num_pad = num_segment* M -length(signal);
paddedSignal = padarray(signal,[num_pad 0],0,'post');
s = reshape(paddedSignal,M,num_segment)';

% piano_fft = zeros(size(s));
% for ss = 1:num_segment
% piano_fft(ss,:)= fft(s(ss,:));
% end

piano_fft = fft(s');
piano_fft = piano_fft';

error = A.*piano_fft;
error_reshape = reshape(error',[848980 1]);
error_time = ifft(error_reshape')';

piano_hat = P .* piano_fft;
piano_hat_reshape = reshape(piano_hat',[848980 1]);


  sound(abs(error_time))
% sound(abs(piano_hat_reshape))
sound(abs(paddedSignal-error_reshape))
% H_reshape = reshape(H',[848980 1]);
% piano_pad_zero = zeros(numel(a_exp)-length(signal));
% piano_pad = [signal piano_pad_zero];
% piano_fft = fft(piano_pad);
