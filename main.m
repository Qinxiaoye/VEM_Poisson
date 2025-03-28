% Copyright:
% Bingbing Xu
% Contact:
% xubingbingdut@foxmail.com, bingbing.xu@ikm.uni-hannover.de
% Institute of Continuum Mechanics, Leibniz Universität Hannover
% This work is mainly from Yue Yu's previous work https://github.com/Terenceyuyue/mVEM


clear;
load('meshdata32.mat');

k = 1; % order

sumNode = size(node,1);

% global stiffness matrix
[GK] = globalK(node,elem);

% right hand vector
ff = zeros(sumNode,1);

% boundary condition

nodeL = find(node(:,1)<0.001); 
nodeR = find(node(:,1)>0.999); 

uh = zeros(sumNode,1);
uh(nodeL) = 0;
uh(nodeR) = 1;
freeDof = setdiff(1:sumNode,[nodeL;nodeR]');
uh(freeDof) = GK(freeDof,freeDof)\(ff(freeDof)-GK(freeDof,[nodeL;nodeR])*uh([nodeL;nodeR]));

% Plot mesh
subplot(1,2,1);
showmesh(node,elem); 
hold on;
plot(node(:,1),node(:,2),'k.', 'MarkerSize', 4);

% Plot numerical solution
subplot(1,2,2);
showsolution(node,elem,uh(1:size(node,1),1));

