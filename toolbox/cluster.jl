function get_distance(points, point)
    dis = 0.0
    for p in points
        dis += norm(p - point)
    end
    return dis
end

function get_closest_point(points, symmetric_points)
    dis = Inf
    closest_point = symmetric_points[1]
    for point in symmetric_points
        new_dis = get_distance(points, point)
        if new_dis < dis
            dis = new_dis
            closest_point = point
        end
    end
    return closest_point
end