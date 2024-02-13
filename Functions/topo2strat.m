% Convert topography to stratigraphy
function Z=topo2strat(Z);

for k=size(Z,1)-1:-1:1

    Z(k,:)=min(Z(k,:),Z(k+1,:));
   
end
