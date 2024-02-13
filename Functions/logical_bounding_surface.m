function Li_bs = logical_bounding_surface(Li,nt)
%Li = logical array 1 = preserved, 0 = amalgamated

Li_bs = false(size(Li));

for j = 1:size(Li,2)
    %this should find the *endpoints* of continously deposited strata
    Li_bs(1:nt,j) = bwmorph(squeeze(~Li(:,j)),'endpoints');
end
