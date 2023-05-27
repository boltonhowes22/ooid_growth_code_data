function points_fin = normalizePointCloud(x0, y0, z0)
    
    % centroid
    centr = mean([x0, y0, z0]);
    
    % centered points: moves all points so the centroid is [0,0,0]
    points = [x0,y0,z0]-centr;

    % find furthest distance from the centroid
    %furthest_distance = max(sqrt(sum(abs(points).^2)));	
    %furthest_distance = min(sqrt(sum(abs(points).^2)));
    D = max(vecnorm(points - [centr], 2, 2));
    %Ds =  vecnorm(points - [centr], 2, 2);
    
    % scale to furthest point
    %points_fin = points./furthest_distance;
    points_fin = points./D;
    
    %%
   
    % other solution
    
    %bla = 100.*randn(1,10);
    %norm_data = (bla - min(bla)) / ( max(bla) - min(bla) );
    
end
