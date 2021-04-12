clear;clc;
%%Gerar a BER para códigos polares
N=512; K=256; EbN0dBrange=0:1:5; designSNRdB=0;
Ec = 1;
N0 = 2;
EbN0dB = EbN0dBrange;
mu=2;
alpha=2;

MCsize = 2000; %Montecarlo size

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
i = 0;
for j=1:length(EbN0dB)
tt=tic();
    N0 = PCparams.N0;
    Ec = (K/N)*N0*10^(EbN0dB(j)/10);
    
    PCparams.Ec = Ec; %normalized Ec %necessary for pencode(), pdecode()
    %l = 0;
    for l = 1:MCsize % (BER(j)<100)||(BERray(j)<100)
        %%polar + AWGN%%
        u=randi(2,K,1)-1; %Bernoulli(0.5)
        x=pencode(u);
        txvec = (2*x-1)*sqrt(Ec);
        noise = randn(N,1);
%         y = txvec + sqrt(N0/2)*noise;
%         uhat = pdecode(y);
%         
%         nfails = sum(uhat ~= u);
%         
%         BER(j) = BER(j) + nfails;
        %%Polar + AWGN%%%%%
        
        %% Polar + desv + AWGN%%%
        %desv generalizado
        %kappa = 1;
        %eta = 1; 
        %r = alphamu_func(N,alpha,mu); %desvanecimento alpha mu
        %r = kappa_mu_func(N,kappa,mu); %desvanecimento kappa mu
        %r = eta_mu_func(N,eta,mu);    %desvanecimento eta mu
        %fim do desv generalizado
        %y_raychannel = txvec.*r + sqrt(N0/2)*noise; %com desvanecimento
        y_raychannel = txvec + sqrt(N0/2)*noise;
        %y_rayMod = y_raychannel./r;
        uhatRay = pdecode(y_raychannel);
        nfailsray = sum(uhatRay ~= u);
        BERray(j) = BERray(j) + nfailsray;
        %% flags de interrupção
        if BERray(j) > 450 && l > 1600
            break
        end
    end
%     BER(j) = BER(j)/(K*l);
    BERray(j) = BERray(j)/(K*l);

end
 
 h = figure();
 semilogy(EbN0dB,BERray); grid on;
 %titlestr = sprintf('N=%d R=%.2f Polar code performance (designSNR=%.1dB)',N,(K/N),designSNRdB);
 %title("Codificação Polar com AWGN e desvanecimento generalizado");
 legenda = sprintf('ENCOM2020%d bits e alpha = %.2f e mu = %.2f puro',N,alpha,mu);
 %legend("Polar + AWGN + desv");
 legend(legenda);
 xlabel('Eb/N0 (dB)');
 ylabel('Bit Error Rate');
 save(sprintf('%s.mat',legenda),'BERray');
 saveas(h,sprintf('%s.fig',legenda));
