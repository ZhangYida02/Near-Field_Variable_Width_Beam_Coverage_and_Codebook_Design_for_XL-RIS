function angle_matrix = compute_angles(x1, y1, z1, x2, y2, z2, v)

    vector_diff_x = x2 - x1'; 
    vector_diff_y = y2 - y1'; 
    vector_diff_z = z2 - z1'; 

    vector_diff = cat(3, vector_diff_x, vector_diff_y, vector_diff_z); 

    dot_product = sum(vector_diff .* reshape(v, [1, 1, 3]), 3);

    norms_vector_diff = sqrt(sum(vector_diff.^2, 3));
    norms_v = norm(v);  


    cos_theta = dot_product ./ (norms_vector_diff * norms_v);

    cos_theta = max(min(cos_theta, 1), -1);

    angle_matrix = acosd(cos_theta); 

    angle_matrix = min(angle_matrix, 180 - angle_matrix);
end
