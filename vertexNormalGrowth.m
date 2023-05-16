function [new_V] = vertexNormalGrowth(v, increment)
% INPUT 
% v: verteces
% increment: how much to grow -- should be the max difference between
%            corteces 
%
%
% OUTPUT
% v_new: the new verteces 
%
%
%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%
% EXAMPLE:
% [v, ~, ~, ~] = stlReadAscii('/Users/boltonhowes/Dropbox (Princeton)/ooids/Ooid3DModels/ooids_gatesburg/Ooid 1/STL/Mesh gatesburg_ooid1a (Convex hull).stl');
% increment = 5;
% vertexNormalGrowth(v, increment);

% Create Triangulation
DT = delaunayTriangulation(v(:,1),v(:,2),v(:,3)); %Note: Shape will always be a convex triangulation
[T,Xb] = freeBoundary(DT);
freeBoundaryT = triangulation(T,Xb);

% Find Unit Normals
V = vertexNormal(freeBoundaryT);

% Calculate the new vertex point
new_V = dCalc(V, Xb, increment);


end

function new_p = dCalc(F, p, d)
    
    %a = sqrt(F(1).^2 + F(2).^2 + F(3).^2);
    %b = d/a;
    %b = b*F;
    %new_p = p+b;
    for t = 1:size(F,1)
       
        new_p(t,:) = p(t,:) + (d/sqrt((F(t,1)^2+F(t,2)^2+F(t,3)^2)))*(F(t,:));
    end
end




%%




