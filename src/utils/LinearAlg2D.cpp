//Copyright (c) 2018 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#include "linearAlg2D.h"

#include <cmath> // atan2
#include <cassert>
#include <algorithm> // swap

#include "IntPoint.h" // dot

namespace cura 
{

float LinearAlg2D::getAngleLeft(const Point& a, const Point& b, const Point& c)
{
    Point ba = a - b;
    Point bc = c - b;
    int64_t dott = dot(ba, bc); // dot product
    int64_t det = ba.X * bc.Y - ba.Y * bc.X; // determinant
    float angle = -atan2(det, dott); // from -pi to pi
    if (angle >= 0 )
    {
        return angle;
    }
    else 
    {
        return M_PI * 2 + angle;
    }
}


bool LinearAlg2D::getPointOnLineWithDist(const Point p, const Point a, const Point b, int64_t dist, Point& result)
{
    //         result
    //         v
    //   b<----r---a.......x
    //          '-.        :
    //              '-.    :
    //                  '-.p
    const Point ab = b - a;
    const int64_t ab_size = vSize(ab);
    const Point ap = p - a;
    const int64_t ax_size = (ab_size < 50)? dot(normal(ab, 1000), ap) / 1000 : dot(ab, ap) / ab_size;
    const int64_t ap_size2 = vSize2(ap);
    const int64_t px_size = sqrt(std::max(int64_t(0), ap_size2 - ax_size * ax_size));
    if (px_size > dist)
    {
        return false;
    }
    const int64_t xr_size = sqrt(dist * dist - px_size * px_size);
    if (ax_size <= 0)
    { // x lies before ab
        const int64_t ar_size = xr_size + ax_size;
        if (ar_size < 0 || ar_size > ab_size)
        { // r lies outisde of ab
            return false;
        }
        else
        {
            result = a + normal(ab, ar_size);
            return true;
        }
    }
    else if (ax_size >= ab_size)
    { // x lies after ab
        //         result
        //         v
        //   a-----r-->b.......x
        //          '-.        :
        //              '-.    :
        //                  '-.p
        const int64_t ar_size = ax_size - xr_size;
        if (ar_size < 0 || ar_size > ab_size)
        { // r lies outisde of ab
            return false;
        }
        else
        {
            result = a + normal(ab, ar_size);
            return true;
        }
    }
    else // ax_size > 0 && ax_size < ab_size
    { // x lies on ab
        //            result is either or
        //         v                       v
        //   a-----r-----------x-----------r----->b
        //          '-.        :        .-'
        //              '-.    :    .-'
        //                  '-.p.-'
        //           or there is not result:
        //         v                       v
        //         r   a-------x---->b     r
        //          '-.        :        .-'
        //              '-.    :    .-'
        //                  '-.p.-'
        // try r in both directions
        const int64_t ar1_size = ax_size - xr_size;
        if (ar1_size >= 0)
        {
            result = a + normal(ab, ar1_size);
            return true;
        }
        const int64_t ar2_size = ax_size + xr_size;
        if (ar2_size < ab_size)
        {
            result = a + normal(ab, ar2_size);
            return true;
        }
        return false;
    }
}

bool LinearAlg2D::lineSegmentsCollide(Point a_from_transformed, Point a_to_transformed, Point b_from_transformed, Point b_to_transformed)
{
    assert(std::abs(a_from_transformed.Y - a_to_transformed.Y) < 2 && "line a is supposed to be transformed to be aligned with the X axis!");
    assert(a_from_transformed.X - 2 <= a_to_transformed.X && "line a is supposed to be aligned with X axis in positive direction!");
    if ((b_from_transformed.Y >= a_from_transformed.Y && b_to_transformed.Y <= a_from_transformed.Y) || (b_to_transformed.Y >= a_from_transformed.Y && b_from_transformed.Y <= a_from_transformed.Y))
    {
        if(b_to_transformed.Y == b_from_transformed.Y)
        {
            if (b_to_transformed.X < b_from_transformed.X)
            {
                std::swap(b_to_transformed.X, b_from_transformed.X);
            }
            if (b_from_transformed.X > a_to_transformed.X)
            {
                return false;
            }
            if (b_to_transformed.X < a_from_transformed.X)
            {
                return false;
            }
            return true;
        }
        else
        {
            int64_t x = b_from_transformed.X + (b_to_transformed.X - b_from_transformed.X) * (a_from_transformed.Y - b_from_transformed.Y) / (b_to_transformed.Y - b_from_transformed.Y);
            if (x >= a_from_transformed.X && x <= a_to_transformed.X)
            {
                return true;
            }
        }
    }
    return false;
}

Point LinearAlg2D::intersection(LineSegment a, LineSegment b)
{
    auto det = [](Point a, Point b) // determinant
        {
            return a.X * b.Y - a.Y * b.X;
        };
    Point a_vec = a.getVector();
    Point b_vec = b.getVector();
    coord_t qwe = det(b.from - a.from, a_vec); // TODO: understand formula and use understandable variable names
    coord_t asd = det(a_vec, b_vec);
    if (asd == 0)
    { // lineas are parallel or collinear
        return (a.from + a.to + b.from + b.to) / 4; // TODO: this is inaccurate
    }
    return b.from + b_vec * qwe / asd;
}

bool LinearAlg2D::areParallel(LineSegment a, LineSegment b, coord_t allowed_error)
{
    Point a_vec = a.getVector();
    Point b_vec = b.getVector();
    coord_t a_size = vSize(a_vec);
    coord_t b_size = vSize(b_vec);
    if (a_size == 0 || b_size == 0)
    {
        return true;
    }
    coord_t dot_size = std::abs(dot(a_vec, b_vec));
    coord_t dot_diff = std::abs(dot_size - a_size * b_size);
    coord_t allowed_dot_error = allowed_error * std::sqrt(dot_size);
    return dot_diff < allowed_dot_error;
}

bool LinearAlg2D::areCollinear(LineSegment a, LineSegment b, coord_t allowed_error)
{
    bool lines_are_parallel = areParallel(a, b, allowed_error);
    bool to_b_from_is_on_line = areParallel(LineSegment(a.from, b.from), b, allowed_error);
    bool to_b_to_is_on_line = areParallel(LineSegment(a.from, b.to), b, allowed_error);
    return lines_are_parallel && to_b_from_is_on_line && to_b_to_is_on_line;
}

coord_t LinearAlg2D::projectedLength(LineSegment to_project, LineSegment onto)
{
    const Point& a = to_project.from;
    const Point& b = to_project.to;
    const Point& c = onto.from;
    const Point& d = onto.to;
    Point cd = d - c;
    coord_t cd_size = vSize(cd);
    assert(cd_size > 0);
    Point ca = a - c;
    coord_t a_projected = dot(ca, cd);
    Point cb = b - c;
    coord_t b_projected = dot(cb, cd);
    return (b_projected - a_projected) / cd_size;
}

LineSegment LinearAlg2D::project(LineSegment to_project, LineSegment onto)
{
    return LineSegment(project(to_project.from, onto), project(to_project.to, onto));
}

Point LinearAlg2D::project(Point p, LineSegment onto)
{
    const Point& a = onto.from;
    const Point& b = onto.to;
    Point ab = b - a;
    coord_t ab_size = vSize(ab);
    assert(ab_size > 0);
    Point pa = p - onto.from;
    coord_t projected_length = dot(ab, pa) / ab_size;
    return onto.from + normal(ab, projected_length);
}

coord_t LinearAlg2D::getTriangleArea(Point a, Point b)
{
    return std::abs(a.X * b.Y - a.Y * b.X) / 2;
}

} // namespace cura
