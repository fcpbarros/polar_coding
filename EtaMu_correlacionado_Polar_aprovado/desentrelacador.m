function [output] = desentrelacador(data,key,depth)

length_orig = length(data);
blocks = ceil(length_orig/depth);
%% Adicionar zeros para deixar compat�vel com qualquer profundidade
%% Os zeros s�o adicionados no final
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
invKey(key) = 1:blocks;
%% Entrela�ar
auxiliar = copia_reshape(:,invKey);
output = reshape(auxiliar,size(copia));
output(isnan(output)) = [];

end

