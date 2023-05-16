function [L2err,lmcosi, diff_cosi, tracker] = fitAllBands(verts, L)
lmcosi = [];
tracker = [];
diff_cosi = [];
for this_ooid = 1:length(verts) 
    for this_band = 1:length(verts{this_ooid})
        
        x = verts{this_ooid}{this_band}(:,1);
        y = verts{this_ooid}{this_band}(:,2);
        z = verts{this_ooid}{this_band}(:,3);

        if this_ooid == 1 && this_band ==1 
            [L2err(this_ooid, this_band), lmcosi{this_ooid}{this_band}]=fitOoid([x, y, z], L);
        else
            [L2err(this_ooid, this_band), lmcosi{this_ooid}{this_band}]=fitOoid([x, y, z], L);
        end
        
        %{
        if this_band ~= 1
            % so we want to find the difference between the 1 and 2 bands
         
            diff_cosi(:,:,end+1) = lmcosi(:,:,end-1) - lmcosi(:,:,end);
                        
        end
        %}
        

        tracker(end+1, this_ooid, this_band) = this_ooid;
        tracker(end+1, this_ooid, this_band) = this_band;
    end
end
end