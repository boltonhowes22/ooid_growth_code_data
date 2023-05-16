function npc = normalizeLoop(vert_struct)
% Function to normalize all point clouds. 
% Makes this size-invariant 


    for this_ooid = 1:length(vert_struct)
        for this_band = 1:length(vert_struct{this_ooid})
            verts = vert_struct{this_ooid}{this_band};

            npc{this_ooid}{this_band} = normalizePointCloud(verts(:,1),...
                verts(:,2), verts(:,3));
            
        end
    end

end
