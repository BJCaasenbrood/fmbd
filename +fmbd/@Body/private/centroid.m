function c = centroid(X, T)
% Volume centroid (straight-sided tets).
% Reference: https://en.wikipedia.org/wiki/Tetrahedron#Centroid
% 
a = X(T(:,2),:) - X(T(:,1),:);  
b = X(T(:,3),:) - X(T(:,1),:); 
v = dot(cross(a, b, 2), X(T(:,4),:) - X(T(:,1),:), 2);
c = v'*(X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:) + X(T(:,4),:))/4/sum(v);
end
