clear;clc;
data = rand(1,1024);
depth = 8;

[data_entre,perm] = entrelacador(data,depth);

data_desentre = desentrelacador(data_entre,perm,depth);

data_desentre - data