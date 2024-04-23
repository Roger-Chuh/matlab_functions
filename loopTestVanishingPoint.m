function loopTestVanishingPoint()

fail_use_vec3_err = 0;
fail_use_dist_err = 0;

for i = 1 : 50
    try
        testVanishingPoint(true);
    catch
        fail_use_vec3_err = fail_use_vec3_err + 1;
    end
end


for i = 1 : 50
    try
        testVanishingPoint(false);
    catch
        fail_use_dist_err = fail_use_dist_err + 1;
    end
end

fail_use_vec3_err
fail_use_dist_err

end