data = [0 1 0 0 1 0 1 0 0 1 1 1 1 1 0 0];
depth = 3;
length_orig = length(data);
blocks = ceil(length_orig/depth);
%% Adicionar NaNs para deixar compatível com qualquer profundidade
%% Os NaNs são adicionados no final
copia = data;
if mod(length_orig,depth) > 0
   nan_to_add = depth - mod(length_orig,depth);
   add_nan = nan(1,nan_to_add);
   copia = [data add_nan];
end
%% Reformulação para blocos de tamanho [blocks depth]
%% As palavras estarão nas colunas
copia_reshape = reshape(copia,[depth blocks]);
%% Gerar os coeficientes de permutação
%perm = randperm(blocks);
perm = [4,5,6,3,1,2];
%% Entrelaçar
auxiliar = copia_reshape(:,perm);
output = reshape(auxiliar,size(copia));
output(isnan(output)) = [];