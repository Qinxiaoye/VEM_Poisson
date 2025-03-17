clear;
pp = load("NLIST.DAT");
pp(:,1) = [];
tt = load("ELIST.DAT");
tt(:,1:1) = [];

[node,elem] = dualMesh(pp(:,1:2),tt(:,1:3));

showmesh(node,elem);