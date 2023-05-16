function rot_inv = rotInv(lmcosi)
% Rotation-invariant shape feature F(l)
% Input just the cos, and sin terms 

    %tracker_begin = [1,2,4,7,11,16];
    %tracker_end = [1,3,6,10,15,21];
    
    tracker_end = 1;
    tracker_begin = 1;
    x = 2;
    while tracker_end(end) < size(lmcosi,1)
       
        tracker_end(x) = x + tracker_end(x-1);
        tracker_begin(x) = x-1 + tracker_begin(x-1);
        x =x +1;
    end
    
    
   
    

    
    for this_degree = 1:length(tracker_begin)
        
        %t = lmcosi(tracker_begin(this_degree):tracker_end(this_degree),1:2).^2;
        t = abs(lmcosi(tracker_begin(this_degree):tracker_end(this_degree),1:2));
        rot_inv(this_degree) = sum(t(:));
        
        if sum(t(:)) < 0 
            keyboard
        end
    end
    
end