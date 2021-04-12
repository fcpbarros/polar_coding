function [output] = desentrelacador(data,key,depth)

length_orig = length(data);
blocks = ceil(length_orig/depth);
%% Adicionar zeros para deixar compatível com qualquer profundidade
%% Os zeros são adicionados no final
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
invKey(key) = 1:blocks;
%% Entrelaçar
auxiliar = copia_reshape(:,invKey);
output = reshape(auxiliar,size(copia));
output(isnan(output)) = [];

end

