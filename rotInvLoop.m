function rot_inv = rotInvLoop(lmcosi)

    for this_ooid = 1:length(lmcosi)
        for this_band = 1:length(lmcosi{this_ooid})
            
            rot_inv{this_ooid}{this_band}= rotInv(lmcosi{this_ooid}{this_band}(:, 3:4));
            
        end
    end

end