clear;clc;
%%Gerar a BER para códigos polares
N=1024; K=512; EbN0dBrange=0:1:8; designSNRdB=0;
Ec = 1;
N0 = 2;
EbN0dB = EbN0dBrange;
mu=1;                                       %mu = 1; alpha = 2 -> Rayleigh
alpha = 2;
fd = 5;
fs = 60; % Sampling frequency [Hz]
Var = 1; % Variance of the Eta-Mu signal / 0 < Var < Infinity
Lambda = 0.8;

D=N; %Interleaver Depth. More the value of D more is the randomness

MCsize = 1000; %Montecarlo size

global PCparams;

PCparams = struct('N', N, ...
    'K', K, ...
    'n', 0, ...
    'FZlookup', zeros(N,1), ...
    'Ec', Ec, ...
    'N0', N0, ...        %N0/2 =1 normalization
    'LLR', zeros(1,2*N-1), ...
    'BITS', zeros(2,N-1), ...
    'designSNRdB', designSNRdB);

PCparams.n = log2(N);

PCparams.FZlookup = pcc(N,K,designSNRdB);

BER = zeros(1,length(EbN0dB));
BERray = zeros(1,length(EbN0dB));

for j=1:length(EbN0dB)
    tt=tic();
    N0 = PCparams.N0;
    Ec = (K/N)*N0*10^(EbN0dB(j)/10);
    
    PCparams.Ec = Ec; %normalized Ec %necessary for pencode(), pdecode()
    %l = 0;
    for l = 1:MCsize % (BER(j)<100)||(BERray(j)<100)
        %%polar + AWGN%%
        u=randi(2,1,K)-1; %Bernoulli(0.5)
        x=pencode(u);
        %% entrelaçador
        %Make the length of the data to be a multiple of D. This is for
        %demonstration only.Otherwise the empty spaces has to be zero filled.
        x=repmat(x,1,1); %send n blocks of specified data pattern
        if mod(length(x),D) ~=0
            x=[x x(1:(D*(fix(length(x)/D+1))-length(x)))];
        end
        
        %Writing into the interleaver row-by-row
        permuterIndex=randperm(D);
        intrlvrOutput=[];
        index=1;
        
        for i=1:fix(length(x)/D)
            intrlvrOutput=[intrlvrOutput x(permuterIndex+(i-1)*D)];
        end
        
        for i=1:mod(length(x),D)
            intrlvrOutput=[intrlvrOutput x(permuterIndex(i)+fix(length(x)/D)*D)];
        end
        
        intrlvrOutput = reshape(intrlvrOutput,[N,1]);
        
        %x = randintrlv(x,8); %%entrelaçador
        %% fim do entrelaçador
        %%
        txvec = (2*intrlvrOutput-1)*sqrt(Ec);
        noise = randn(N,1);
        r = alphamu_corre_coefs(N,mu,alpha,Var,Lambda,fd,fs);    %desvanecimento eta mu
        %fim do desv generalizado
        y = txvec.*r + sqrt(N0/2)*noise;
        %y_rayMod = y_raychannel./r;
        %% desentrelaçador
        
        %y = reshape(y,[D,N/D]);
        
        deintrlvrOutput=[];
        
        for i=1:fix(length(y)/D)
            deintrlvrOutput(permuterIndex+(i-1)*D)=y((i-1)*D+1:i*D);
        end
        
        for i=1:mod(length(y),D)
            deintrlvrOutput((fix(length(y)/D))*permuterIndex+i)=y((i+1):(i+1)*D);
        end
        
        
        %y = randdeintrlv(y,8);
        
        
        %% fim do desentrelaçador
        uhatRay = pdecode(deintrlvrOutput);
        nfailsray = sum(uhatRay' ~= u);
        BERray(j) = BERray(j) + nfailsray;
        %% flags de interrupção
        if BERray(j) > 350 && l > 650
            break
        end
    end
    %     BER(j) = BER(j)/(K*l);
    BERray(j) = BERray(j)/(K*l);
    
end

h = figure();
semilogy(EbN0dB,BERray); grid on;
%titlestr = sprintf('N=%d R=%.2f Polar code performance (designSNR=%.1dB)',N,(K/N),designSNRdB);
%title("Codificação Polar com AWGN e desvanecimento generalizado correl entrelacado");
legenda = sprintf('ENTRELACADOR_NOVO%d bits e alpha = %.2f e mu = %.2f 8corr10_entr_%d',N,alpha,mu,D);
%legend("Polar + AWGN + desv");
legend(legenda);
xlabel('Eb/N0 in dB');
ylabel('Bit Error Rate');
save(sprintf('%s.mat',legenda),'BERray')
saveas(h,sprintf('%s.fig',legenda));
