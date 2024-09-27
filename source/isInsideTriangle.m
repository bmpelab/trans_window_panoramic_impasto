function inside = isInsideTriangle(point, A, B, C)
    v0 = C - A;
    v1 = B - A;
    v2 = point - A;

    dot00 = dot(v0, v0);
    dot01 = dot(v0, v1);
    dot02 = dot(v0, v2);
    dot11 = dot(v1, v1);
    dot12 = dot(v1, v2);

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    inside = (u >= 0) & (v >= 0) & (u + v <= 1);
end