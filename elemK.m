function AK = elemK(elemID,centroid,diameter,node,elem)
% calculate elem stiffness matrix
% for k = 1
% input: elemID,centroid,diameter,node,elem
% output: AK

k = 1;
index = elem{elemID};
Nv = length(index);
hK = diameter(elemID);

x = node(index ,1); y = node(index ,2);

% calculate D
D = myM([x,y],k,hK,centroid(elemID,:)); % Eq.(58)

% calculate B
Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]]; % k = 1
rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
normVec = 0.5*[y(rotid2)-y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
B = Gradm*normVec; % B, Eq.(60)

% constraint
Bs = B;  Bs(1,:) = 1/Nv; % constraint, Eq.(37)
% consistency relation
G = B*D;  Gs = Bs*D;  % Eq.(61)

% --------- local stiffness matrix ---------
Pis = Gs\Bs;  
Pi  = D*Pis;
I = eye(size(Pi));
AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);