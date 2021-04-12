data = [0 1 0 0 1 0 1 0 0 1 1 1 1 1 0 0];
depth = 3;
length_orig = length(data);
blocks = ceil(length_orig/depth);
%% Adicionar NaNs para deixar compat�vel com qualquer profundidade
%% Os NaNs s�o adicionados no final
copia = data;
if mod(length_orig,depth) > 0
   nan_to_add = depth - mod(length_orig,depth);
   add_nan = nan(1,nan_to_add);
   copia = [data add_nan];
end
%% Reformula��o para blocos de tamanho [blocks depth]
%% As palavras estar�o nas colunas
copia_reshape = reshape(copia,[depth blocks]);
%% Gerar os coeficientes de permuta��o
%perm = randperm(blocks);
perm = [4,5,6,3,1,2];
%% Entrela�ar
auxiliar = copia_reshape(:,perm);
output = reshape(auxiliar,size(copia));
output(isnan(output)) = [];