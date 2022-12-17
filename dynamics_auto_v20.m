function X_dot = dynamics_auto_v20(K,A,t,X,K_id_total,p_id_total,n_id_total)

%make p1
p1 = ones(size(K,1),1);
box = find(p_id_total(:,1));
p1(box) = X(nonzeros(p_id_total(:,1)));

%make p2
p2 = ones(size(K,1),1);
box = find(p_id_total(:,2));
p2(box) = X(nonzeros(p_id_total(:,2)));

%make p3
p3 = ones(size(K,1),1);
box = find(p_id_total(:,3));
p3(box) = X(nonzeros(p_id_total(:,3)));

%make p4
p4 = ones(size(K,1),1);
box = find(p_id_total(:,4));
p4(box) = X(nonzeros(p_id_total(:,4)));

%make n1
n1 = ones(size(K,1),1);
box = find(n_id_total(:,1));
n1(box) = X(nonzeros(n_id_total(:,1)));

%make n2
n2 = ones(size(K,1),1);
box = find(n_id_total(:,2));
n2(box) = X(nonzeros(n_id_total(:,2)));

%make n3
n3 = ones(size(K,1),1);
box = find(n_id_total(:,3));
n3(box) = X(nonzeros(n_id_total(:,3)));

%make n4
n4 = ones(size(K,1),1);
box = find(n_id_total(:,4));
n4(box) = X(nonzeros(n_id_total(:,4)));

V = K_id_total(:,1).* K(:,1) .* p1 .*p2.* p3.* p4 + K_id_total(:,2) .* n1 .* n2 .* n3 .* n4;
X_dot =A*V;

end