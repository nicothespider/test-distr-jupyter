function x=signalgen(type, Wn, N, amplitude, options)
% X = signalgen(TYPE, Wn, N, LEVELS)

if nargin<5
    options = [];
end
if nargin<4 || isempty(amplitude)
    amplitude = [-1 1];
end
if nargin<3  || isempty(N)
    N = 200;
end
if nargin<2  || isempty(Wn)
    Wn = [0 1];
end
if length(Wn)==1
    Wn = [0 Wn];
end
N = ceil(N);
freqmin = Wn(1) / 2; % (1 pour la fréquence de Nyquist)
freqmax = Wn(2) / 2; 

switch(lower(type))
    case {'sga', 'rgs'}
        x = randn(1,N);
        if Wn(1)~=0 || Wn(2)~=1
            if Wn(1)==0
                [B,A] = butter(8,Wn(2));
            elseif Wn(2)==1
                [B,A] = butter(8,Wn(1),'high');
            else
                [B,A] = butter(8,Wn);
            end
        end
        x = filtfilt(B,A,x); 
    case {'sba', 'rbs'}
        x = randn(1,N);
        if Wn(1)~=0 || Wn(2)~=1
            if Wn(1)==0
                [B,A] = butter(8,Wn(2));
            elseif Wn(2)==1
                [B,A] = butter(8,Wn(1),'high');
            else
                [B,A] = butter(8,Wn);
            end
        end
        x = filtfilt(B,A,x);
        x = sign(x);
    case {'sin', 'sine'}
        % options = [nsin]
        if isempty(options)
            nsin = 30;
        else
            nsin = options(1);
        end
        % Génération
        t = 1:1:N;
        x = zeros(size(t));
        for i=1:nsin
            f = freqmin+(freqmax-freqmin)*i/nsin;
            phi = 2*pi*rand(1);
            x = x + sin(2*pi*f*t + phi);
        end
    case {'sbpa', 'prbs'}
        % Positions
        feedback = hex2dec([...
            '0000000000',
            '0000000003',
            '0000000005',
            '0000000009',
            '0000000012',
            '0000000021',
            '0000000041',
            '000000008e',
            %'00000000b1',
            '0000000108',
            '0000000204',
            '0000000402',
            '0000000829',
            '000000100d',
            '0000002015',
            '0000004001',
            '0000008016',
            '0000010004',
            '0000020013',
            '0000040013',
            '0000080004',
            '0000100002',
            '0000200001',
            '0000400010',
            '000080000d',
            '0001000004',
            '0002000023',
            '0004000013',
            '0008000004',
            '0010000002',
            '0020000029',
            '0040000004',
            '0080000057',
            '0100000029',
            '0200000073',
            '0400000002',
            '080000003b',
            '100000001f',
            '2000000031',
            '4000000008'
            ]); 
        T = floor(1/freqmax);
        %n = ceil(log2(N/T+1)); % par excès
        %n = floor(log2(N/T+1)); % par défaut
        n = round(log2(N/T+1)); % par arrondi
        n = min([n, length(feedback)]);
        feedback = str2num(char(cellstr(dec2bin(feedback(n),n)')))';
        % Génération
        %%%len = 2^n-1;
        len = ceil(N/T);
        x = ones(1,len+n);
        xx = zeros(1,N);
        for i=1+n:1:len+n
            x(i) = mod(sum(and(feedback,x(i-n:i-1))),2);
            xx((i-n-1)*T+1:(i-n)*T) = ones(1,T)*x(i);
        end
%         x = filter(ones(1,T),1,upsample(x,T));
%         x = x(n+1:N+n);
        x = xx;

end

% Mise à l'échelle
x = amplitude(1) + (amplitude(2)-amplitude(1))*(x-min(x))/(max(x)-min(x));
