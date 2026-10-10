/*
 * Geometry of the cells cut by an interface, independent of the mesh:
 *  - polygons cut once (cells of a 2D mesh, faces of a 3D mesh in the
 *    coordinates of their plane): position of the interface (polyline,
 *    polynomial or implicit graph above its chord) and quadratures of the two
 *    parts by sectors from a point with respect to which they are star-shaped;
 *  - tetrahedra cut by the interface, approximated by subdivision (fallback
 *    of the 3D cut cells).
 * The level set functions given to the polygon algorithms are expressed in
 * the coordinates of the plane (see planar_level_set).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "diskpp/mesh/point.hpp"
#include "diskpp/mesh/cut_types.hpp"
#include "diskpp/mesh/level_set.hpp"
#include "diskpp/quadratures/quadratures.hpp"

namespace disk {
namespace cutgeo {

template<typename T>
using vec2 = Eigen::Matrix<T, 2, 1>;

template<typename T>
using vec3 = Eigen::Matrix<T, 3, 1>;

template<typename T>
vec2<T> to_vec(const point<T,2>& p) { return vec2<T>(p.x(), p.y()); }

template<typename T>
T cross(const point<T,2>& u, const point<T,2>& v) { return u.x()*v.y() - u.y()*v.x(); }

template<typename T>
T cross(const point<T,2>& u, const vec2<T>& v) { return u.x()*v(1) - u.y()*v(0); }

/* exact equality of two points (copies of the same point) */
template<typename T, size_t DIM>
bool
same_point(const point<T,DIM>& a, const point<T,DIM>& b) {
    for (size_t i = 0; i < DIM; i++)
        if (a[i] != b[i])
            return false;
    return true;
}

/* Gauss-Legendre rule with n points on [0,1] (Newton iterations on P_n) */
template<typename T>
const std::vector<std::pair<T, T>>&
gauss_legendre_01(size_t n) {

    static std::map<size_t, std::vector<std::pair<T, T>>> cache;
    auto it = cache.find(n);
    if (it != cache.end())
        return it->second;

    std::vector<std::pair<T, T>> rule(n);
    for (size_t i = 0; i < n; i++) {
        T x = std::cos(std::acos(T(-1)) * (i + 0.75) / (n + 0.5));
        T p0 = 1, p1 = x, dp = 1;
        for (size_t iter = 0; iter < 100; iter++) {
            p0 = 1; p1 = x;
            for (size_t k = 2; k <= n; k++) {
                T p2 = ((2.0*k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
                p0 = p1; p1 = p2;
            }
            dp = n * (x * p1 - p0) / (x * x - 1.0);
            T dx = p1 / dp;
            x -= dx;
            if (std::abs(dx) < 1e-16)
                break;
        }
        p0 = 1; p1 = x;
        for (size_t k = 2; k <= n; k++) {
            T p2 = ((2.0*k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
            p0 = p1; p1 = p2;
        }
        dp = (n == 1) ? 1.0 : n * (x * p1 - p0) / (x * x - 1.0);
        rule[i] = std::make_pair(0.5 * (1.0 - x), 1.0 / ((1.0 - x * x) * dp * dp));
    }
    cache[n] = rule;
    return cache[n];
}

/* Rule of degree "degree" on the reference triangle (0,0), (1,0), (0,1)
 * (weights summing to 1/2). The weights of a mapped triangle are scaled by
 * the cross product of its edges: this works for very thin triangles, where
 * the Kahan formula of the area used by triangle_gauss() returns NaN. */
template<typename T>
const std::vector<std::pair<point<T,2>, T>>&
reference_triangle_rule(size_t degree) {

    static std::map<size_t, std::vector<std::pair<point<T,2>, T>>> cache;
    auto it = cache.find(degree);
    if (it != cache.end())
        return it->second;

    std::vector<std::pair<point<T,2>, T>> rule;
    auto qps = disk::quadrature::triangle_gauss(degree, point<T,2>(0.0, 0.0), point<T,2>(1.0, 0.0),
                                                point<T,2>(0.0, 1.0));
    for (auto& qp : qps)
        rule.push_back(std::make_pair(qp.point(), qp.weight()));
    cache[degree] = rule;
    return cache[degree];
}

/* Collapsed Gauss rule of degree "degree" on the reference tetrahedron
 * (positive weights summing to 1/6) */
template<typename T>
const std::vector<std::pair<std::array<T,3>, T>>&
reference_tetrahedron_rule(size_t degree) {

    static std::map<size_t, std::vector<std::pair<std::array<T,3>, T>>> cache;
    auto it = cache.find(degree);
    if (it != cache.end())
        return it->second;

    /* x = a, y = b (1-a), z = c (1-a)(1-b), Jacobian (1-a)^2 (1-b) */
    size_t n = (degree + 2)/2 + 1;
    const auto& g = gauss_legendre_01<T>(n);
    std::vector<std::pair<std::array<T,3>, T>> rule;
    for (auto& ga : g)
        for (auto& gb : g)
            for (auto& gc : g) {
                T a = ga.first, b = gb.first, c = gc.first;
                std::array<T,3> x = {a, b*(1-a), c*(1-a)*(1-b)};
                rule.push_back(std::make_pair(x, ga.second*gb.second*gc.second*(1-a)*(1-a)*(1-b)));
            }
    cache[degree] = rule;
    return cache[degree];
}

/* Root of f closest to s0 in [smin, smax], f(s) and f'(s) being given by
 * fdf(s) -> (f, f'). The root is bracketed by steps growing geometrically
 * from s0 (the first step is "step"), then refined by Newton iterations
 * safeguarded by bisection. */
template<typename T, typename FDF>
std::optional<T>
find_root_near(const FDF& fdf, T smin, T smax, T s0, T step, T tol)
{
    auto f0 = fdf(s0).first;
    if (f0 == 0)
        return s0;

    bool found = false;
    T a = s0, b = s0, fa = f0, fb = f0;
    T lo = s0, flo = f0, hi = s0, fhi = f0;
    for (size_t k = 0; k < 80 && !found; k++) {
        T d = step * std::pow(T(2), T(k));
        bool more = false;
        /* candidates on both sides at this distance */
        std::optional<std::array<T,4>> right, left;
        if (hi < smax) {
            T s = std::min(s0 + d, smax);
            T fs = fdf(s).first;
            if ((fhi < 0) != (fs < 0) || fs == 0)
                right = std::array<T,4>{hi, s, fhi, fs};
            hi = s; fhi = fs;
            more = true;
        }
        if (lo > smin) {
            T s = std::max(s0 - d, smin);
            T fs = fdf(s).first;
            if ((flo < 0) != (fs < 0) || fs == 0)
                left = std::array<T,4>{s, lo, fs, flo};
            lo = s; flo = fs;
            more = true;
        }
        if (right && left) {
            /* the closest root, estimated by linear interpolation */
            auto est = [&](const std::array<T,4>& br) {
                return std::abs(br[0] - br[2]*(br[1] - br[0])/(br[3] - br[2]) - s0);
            };
            if (est(*left) < est(*right))
                right.reset();
        }
        if (right) {
            a = (*right)[0]; b = (*right)[1]; fa = (*right)[2]; fb = (*right)[3];
            found = true;
        }
        else if (left) {
            a = (*left)[0]; b = (*left)[1]; fa = (*left)[2]; fb = (*left)[3];
            found = true;
        }
        if (!more)
            break;
    }
    if (!found)
        return std::nullopt;
    if (fa == 0)
        return a;
    if (fb == 0)
        return b;

    T s = a - fa*(b - a)/(fb - fa);
    for (size_t it = 0; it < 100; it++) {
        auto fd = fdf(s);
        if (fd.first == 0)
            return s;
        if ((fd.first < 0) == (fa < 0)) {
            a = s; fa = fd.first;
        }
        else {
            b = s; fb = fd.first;
        }
        T sn = s;
        bool newton = (fd.second != 0);
        if (newton) {
            sn = s - fd.first/fd.second;
            newton = (sn > std::min(a, b) && sn < std::max(a, b));
        }
        if (!newton)
            sn = 0.5*(a + b);
        bool done = std::abs(sn - s) <= tol || std::abs(b - a) <= tol;
        s = sn;
        if (done)
            break;
    }
    return s;
}

/* Level set of a DIM-dimensional space restricted to the plane
 * x = origin + U u (2D meshes: identity), in the coordinates u of the plane */
template<typename T, size_t DIM>
class planar_level_set {

    const level_set<T, DIM>&    m_ls;
    point<T, DIM>               m_origin;
    Eigen::Matrix<T, DIM, 2>    m_U;
    bool                        m_identity;

public:
    planar_level_set(const level_set<T, DIM>& ls)
        : m_ls(ls), m_U(Eigen::Matrix<T, DIM, 2>::Identity()), m_identity(true)
    {
        static_assert(DIM == 2, "the identity frame is only available in 2D");
    }

    planar_level_set(const level_set<T, DIM>& ls, const point<T, DIM>& origin, const Eigen::Matrix<T, DIM, 2>& U)
        : m_ls(ls), m_origin(origin), m_U(U), m_identity(false)
    {}

    point<T, DIM>
    to_global(const point<T,2>& u) const {
        if constexpr (DIM == 2) {
            if (m_identity)
                return u;
        }
        Eigen::Matrix<T, DIM, 1> x = m_origin.to_vector() + m_U * to_vec(u);
        return point<T, DIM>(x);
    }

    Eigen::Matrix<T, DIM, 1>
    to_global_vector(const vec2<T>& v) const {
        return m_U * v;
    }

    point<T,2>
    to_local(const point<T, DIM>& x) const {
        if constexpr (DIM == 2) {
            if (m_identity)
                return x;
        }
        vec2<T> u = m_U.transpose() * (x - m_origin).to_vector();
        return point<T,2>(u(0), u(1));
    }

    T operator()(const point<T,2>& u) const { return m_ls(to_global(u)); }

    vec2<T>
    gradient(const point<T,2>& u) const {
        if constexpr (DIM == 2) {
            if (m_identity)
                return m_ls.gradient(u);
        }
        return m_U.transpose() * m_ls.gradient(to_global(u));
    }
};

/***************************************************************************/
/*                             CUT POLYGONS                                 */
/***************************************************************************/

/* Initialize a polygon cut once by the interface. "vertices" are given
 * counterclockwise; edge i joins the vertices i and i+1, and intersection(i)
 * is the point where the interface crosses it if its ends are on different
 * sides. p0 is the point where the interface leaves the negative vertices (in
 * the counterclockwise order), so that the negative side is on the left of
 * p0 -> p1. Returns false if the polygon is not cut exactly once. */
template<typename T, typename Intersection>
bool
init_cut_polygon(cut_polygon<T>& poly, const std::vector<point<T,2>>& vertices,
                 const std::vector<location>& loc, const std::vector<size_t>& ids,
                 const Intersection& intersection)
{
    size_t n = vertices.size();
    poly.vertices = vertices;
    poly.vertex_loc = loc;
    poly.vertex_ids = ids;

    size_t changes = 0;
    for (size_t i = 0; i < n; i++) {
        auto a = loc[i], b = loc[(i+1)%n];
        if (a == b)
            continue;
        changes++;
        if (a == location::IN_NEGATIVE_SIDE) {
            poly.edge0 = i;
            poly.p0 = intersection(i);
        }
        else {
            poly.edge1 = i;
            poly.p1 = intersection(i);
        }
    }
    if (changes != 2)
        return false;

    poly.h = 0;
    for (size_t i = 0; i < n; i++)
        for (size_t j = i+1; j < n; j++)
            poly.h = std::max(poly.h, (vertices[i] - vertices[j]).to_vector().norm());

    poly.convex = true;
    for (size_t i = 0; i < n; i++) {
        auto e1 = vertices[(i+1)%n] - vertices[i];
        auto e2 = vertices[(i+2)%n] - vertices[(i+1)%n];
        if (cross(e1, e2) < -1e-12 * poly.h * poly.h)
            poly.convex = false;
    }

    poly.kind = interface_kind::polyline;
    poly.polyline = {poly.p0, poly.p1};
    poly.nodes.clear();
    poly.heights.clear();
    return true;
}

/* Polygon of the part "where", the interface being replaced by its chord if
 * chord_only (ProtoN's collect_triangulation_points: the vertices of "where"
 * in counterclockwise order from the first vertex, the interface inserted) */
template<typename T>
std::vector<point<T,2>>
side_polygon(const cut_polygon<T>& poly, location where, bool chord_only = false)
{
    std::vector<point<T,2>> ret;
    const auto& vl = poly.vertex_loc;
    const auto& vs = poly.vertices;
    size_t n = vs.size();

    const std::vector<point<T,2>> chord = {poly.p0, poly.p1};
    const auto& itf = chord_only ? chord : poly.polyline;
    auto insert_interface = [&]() {
        if (where == location::IN_NEGATIVE_SIDE)
            ret.insert(ret.end(), itf.begin(), itf.end());
        else
            ret.insert(ret.end(), itf.rbegin(), itf.rend());
    };

    bool case1 = vl.front() == where && vl.back() != where;
    bool case2 = vl.front() != where && vl.back() == where;
    bool case3 = vl.front() != where && vl.back() != where;

    if (case1 || case2 || case3) {
        for (size_t i = 0; i < n; i++)
            if (vl[i] == where)
                ret.push_back(vs[i]);
        insert_interface();
    }
    else {
        size_t i = 0;
        while (i < n && vl[i] == where)
            ret.push_back(vs[i++]);
        insert_interface();
        while (i < n && vl[i] != where)
            i++;
        while (i < n && vl[i] == where)
            ret.push_back(vs[i++]);
    }

    return ret;
}

/* Is the edge a -> b of a side polygon the chord of the interface? */
template<typename T>
bool
is_chord(const cut_polygon<T>& poly, const point<T,2>& a, const point<T,2>& b, location where)
{
    if (where == location::IN_NEGATIVE_SIDE)
        return same_point(a, poly.p0) && same_point(b, poly.p1);
    return same_point(a, poly.p1) && same_point(b, poly.p0);
}

/* Inward half-planes m.(x - a) >= 0 of the edges of a counterclockwise polygon */
template<typename T>
std::vector<std::pair<vec2<T>, point<T,2>>>
inward_halfplanes(const cut_polygon<T>& poly)
{
    std::vector<std::pair<vec2<T>, point<T,2>>> ret;
    size_t n = poly.vertices.size();
    for (size_t i = 0; i < n; i++) {
        const auto& a = poly.vertices[i];
        const auto& b = poly.vertices[(i+1)%n];
        auto e = (b - a).to_vector();
        T len = e.norm();
        if (len == 0)
            continue;
        ret.push_back(std::make_pair(vec2<T>(-e(1)/len, e(0)/len), a));
    }
    return ret;
}

/* Range [smin, smax] of s such that c + s d lies in the polygon (in its
 * convex hull of the half-planes if convex, in its bounding box otherwise) */
template<typename T>
std::pair<T, T>
line_range(const cut_polygon<T>& poly, const std::vector<std::pair<vec2<T>, point<T,2>>>& halfplanes,
           const point<T,2>& c, const vec2<T>& d)
{
    T smin = -std::numeric_limits<T>::max(), smax = std::numeric_limits<T>::max();
    if (poly.convex) {
        for (auto& hp : halfplanes) {
            T md = hp.first.dot(d);
            T mc = hp.first.dot((c - hp.second).to_vector());
            if (std::abs(md) < 1e-14)
                continue;
            T s = -mc/md;
            if (md > 0)
                smin = std::max(smin, s);
            else
                smax = std::min(smax, s);
        }
        return std::make_pair(smin, smax);
    }
    T b[4] = {poly.vertices[0].x(), poly.vertices[0].x(), poly.vertices[0].y(), poly.vertices[0].y()};
    for (auto& v : poly.vertices) {
        b[0] = std::min(b[0], v.x()); b[1] = std::max(b[1], v.x());
        b[2] = std::min(b[2], v.y()); b[3] = std::max(b[3], v.y());
    }
    for (size_t a = 0; a < 2; a++) {
        T ca = (a == 0) ? c.x() : c.y();
        if (std::abs(d(a)) < 1e-14)
            continue;
        T s1 = (b[2*a] - ca)/d(a), s2 = (b[2*a + 1] - ca)/d(a);
        smin = std::max(smin, std::min(s1, s2));
        smax = std::min(smax, std::max(s1, s2));
    }
    return std::make_pair(smin, smax);
}

/* Is the point in the polygon (with a tolerance)? Convex polygons only;
 * always true for the others. */
template<typename T>
bool
inside(const cut_polygon<T>& poly, const std::vector<std::pair<vec2<T>, point<T,2>>>& halfplanes,
       const point<T,2>& x, T tol)
{
    if (!poly.convex)
        return true;
    for (auto& hp : halfplanes)
        if (hp.first.dot((x - hp.second).to_vector()) < -tol)
            return false;
    return true;
}

/* The interface of a cut polygon as a curve x(t), t in [0,1], from p0 to p1
 * (polynomial and implicit interfaces) */
template<typename T, typename LS2>
class polygon_curve {

    const cut_polygon<T>&                       m_poly;
    const LS2&                                  m_phi;
    vec2<T>                                     m_d, m_tau, m_nu;
    T                                           m_L;
    std::vector<std::pair<vec2<T>, point<T,2>>> m_halfplanes;

public:
    polygon_curve(const cut_polygon<T>& poly, const LS2& phi)
        : m_poly(poly), m_phi(phi)
    {
        m_d   = (poly.p1 - poly.p0).to_vector();
        m_L   = m_d.norm();
        m_tau = m_d / m_L;
        m_nu  = vec2<T>(-m_tau(1), m_tau(0));
        m_halfplanes = inward_halfplanes(poly);
    }

    T length() const { return m_L; }
    const vec2<T>& normal() const { return m_nu; } // left normal of the chord

    point<T,2> chord(T t) const { return m_poly.p0 + (m_poly.p1 - m_poly.p0) * t; }

    /* height of the polynomial interface and its derivative */
    std::pair<T, T>
    height(T t) const
    {
        T h = 0, dh = 0;
        const auto& nodes = m_poly.nodes;
        size_t n = nodes.size();
        for (size_t j = 0; j < n; j++) {
            T lj = 1, dlj = 0;
            for (size_t m = 0; m < n; m++) {
                if (m == j)
                    continue;
                T den = nodes[j] - nodes[m];
                dlj = (dlj * (t - nodes[m]) + lj) / den;
                lj  = lj * (t - nodes[m]) / den;
            }
            h  += m_poly.heights[j] * lj;
            dh += m_poly.heights[j] * dlj;
        }
        return std::make_pair(h, dh);
    }

    /* Root of the level set along nu at the point t of the chord, closest to
     * the chord (implicit interface) */
    std::optional<T>
    implicit_height(T t) const
    {
        auto c = chord(t);
        auto range = line_range(m_poly, m_halfplanes, c, m_nu);
        /* the ends of the chord lie on the boundary: the range is widened */
        T smin = std::min(range.first, T(0)) - T(1e-8)*m_poly.h;
        T smax = std::max(range.second, T(0)) + T(1e-8)*m_poly.h;
        auto fdf = [&](T s) {
            point<T,2> x = c + point<T,2>(m_nu(0)*s, m_nu(1)*s);
            return std::make_pair(m_phi(x), m_phi.gradient(x).dot(m_nu));
        };
        return find_root_near(fdf, smin, smax, T(0), T(1e-3)*m_poly.h, T(1e-15)*m_poly.h);
    }

    /* x(t) and x'(t). Returns false if the implicit interface is not found
     * or is not a graph above the chord at t. */
    bool
    eval(T t, point<T,2>& x, vec2<T>& dx) const
    {
        if (m_poly.kind == interface_kind::polynomial) {
            auto h = height(t);
            x  = chord(t) + point<T,2>(m_nu(0)*h.first, m_nu(1)*h.first);
            dx = m_L * m_tau + h.second * m_nu;
            return true;
        }
        auto s = implicit_height(t);
        if (!s)
            return false;
        x = chord(t) + point<T,2>(m_nu(0) * (*s), m_nu(1) * (*s));
        vec2<T> g = m_phi.gradient(x);
        T gn = g.dot(m_nu);
        if (std::abs(gn) <= 0.2 * g.norm())
            return false;
        dx = m_d - (g.dot(m_d)/gn) * m_nu;
        return true;
    }

    /* Number of Gauss points in t for a rule of degree "degree" along the curve */
    size_t
    nodes_count(size_t degree, bool sector) const
    {
        if (m_poly.kind == interface_kind::polynomial) {
            size_t q = m_poly.nodes.size() - 1;
            return sector ? (degree * q + 2 * q) / 2 + 1 : (degree * q + 2 * q) / 2 + 2;
        }
        return (degree + 1) / 2 + 6;
    }
};

/* Area centroid of a polygon (0 area: vertex average) */
template<typename T>
point<T,2>
polygon_centroid(const std::vector<point<T,2>>& tp, T* area_out = nullptr)
{
    point<T,2> cg;
    T area = 0;
    for (size_t k = 1; k + 1 < tp.size(); k++) {
        auto u = tp[k] - tp[0], v = tp[k+1] - tp[0];
        T a = 0.5 * cross(u, v);
        cg = cg + (tp[0] + tp[k] + tp[k+1]) * (a / 3.0);
        area += a;
    }
    if (area_out)
        *area_out = area;
    if (area == 0) {
        point<T,2> m;
        for (auto& p : tp)
            m = m + p;
        return m / T(tp.size());
    }
    return cg / area;
}

/* Is the part "where" of a cut polygon star-shaped with respect to c? The
 * straight edges and the curve (at nsamples points) are checked. */
template<typename T, typename LS2>
bool
star_shaped(const cut_polygon<T>& poly, const polygon_curve<T, LS2>& crv, location where,
            const point<T,2>& c, size_t nsamples)
{
    T tol = -1e-12 * poly.h * poly.h;
    bool neg = (where == location::IN_NEGATIVE_SIDE);
    auto chord_poly = side_polygon(poly, where, true);
    for (size_t k = 0; k < chord_poly.size(); k++) {
        const auto& a = chord_poly[k];
        const auto& b = chord_poly[(k+1) % chord_poly.size()];
        if (is_chord(poly, a, b, where))
            continue;
        if (cross(a - c, b - c) < tol)
            return false;
    }
    for (size_t k = 0; k <= nsamples; k++) {
        T t = T(k) / T(nsamples);
        point<T,2> x;
        vec2<T> dx;
        if (!crv.eval(t, x, dx))
            return false;
        T det = cross(x - c, dx);  /* negative side: x(t) traversed forward */
        if (!neg)
            det = -det;            /* positive side: backward */
        if (det < tol)
            return false;
    }
    return true;
}

/* Interface of degree "degree" (polynomial interpolation of the heights of
 * the interface above its chord at the Chebyshev-Lobatto nodes, ProtoN's
 * make_curved_interface). The polyline already computed (refine_polyline)
 * is used to check that the curve follows the interface; on success, it is
 * replaced by "samples" points of the curve. Returns false (the polygon
 * keeping its polyline) if the curve is not a graph above the chord in the
 * polygon, leaves it, or if a part of the polygon is not star-shaped with
 * respect to its centroid. */
template<typename T, typename LS2>
bool
make_polynomial_interface(cut_polygon<T>& poly, const LS2& phi, size_t degree, size_t samples = 32)
{
    if (degree < 2)
        return false;

    std::vector<T> nodes(degree + 1);
    for (size_t i = 0; i <= degree; i++)
        nodes[i] = 0.5 * (1.0 - std::cos(std::acos(T(-1)) * T(i) / T(degree)));

    auto saved_kind = poly.kind;
    auto saved_polyline = poly.polyline;
    auto restore = [&]() {
        poly.kind = saved_kind;
        poly.polyline = saved_polyline;
        poly.nodes.clear();
        poly.heights.clear();
        return false;
    };

    T h_cell = poly.h;
    vec2<T> d = (poly.p1 - poly.p0).to_vector();
    T L = d.norm();
    if (!(L > 1e-12 * h_cell))
        return false;
    vec2<T> nu(-d(1) / L, d(0) / L);
    auto halfplanes = inward_halfplanes(poly);

    std::vector<T> heights(degree + 1, 0.0);
    for (size_t i = 1; i < degree; i++) {
        point<T,2> c = poly.p0 + (poly.p1 - poly.p0) * nodes[i];
        auto range = line_range(poly, halfplanes, c, nu);
        T smin = range.first, smax = range.second;
        auto f = [&](T sv) { return phi(c + point<T,2>(nu(0) * sv, nu(1) * sv)); };
        /* bracket the root closest to the chord */
        const size_t N = 64;
        T best_a = 0, best_b = 0, best_dist = std::numeric_limits<T>::max();
        bool found = false;
        T sa = smin, fa = f(sa);
        for (size_t k = 1; k <= N; k++) {
            T sb = smin + (smax - smin) * T(k) / T(N);
            T fb = f(sb);
            if ((fa <= 0 && fb > 0) || (fa > 0 && fb <= 0)) {
                T dist = std::min(std::abs(sa), std::abs(sb));
                if (sa <= 0 && sb >= 0)
                    dist = 0;
                if (dist < best_dist) {
                    best_dist = dist; best_a = sa; best_b = sb; found = true;
                }
            }
            sa = sb; fa = fb;
        }
        if (!found)
            return false;
        T fla = f(best_a);
        for (size_t it = 0; it < 200 && best_b - best_a > 1e-15 * h_cell; it++) {
            T sm = 0.5 * (best_a + best_b);
            T fm = f(sm);
            if ((fla <= 0 && fm > 0) || (fla > 0 && fm <= 0))
                best_b = sm;
            else {
                best_a = sm;
                fla = fm;
            }
        }
        heights[i] = 0.5 * (best_a + best_b);
    }

    poly.kind = interface_kind::polynomial;
    poly.nodes = nodes;
    poly.heights = heights;
    polygon_curve<T, LS2> crv(poly, phi);

    /* the curve must follow the points of the polyline (a kink of the
     * interface, e.g. a corner of the square of Fig. 12, cannot be
     * represented by a polynomial) */
    for (auto& x : saved_polyline) {
        auto r = x - poly.p0;
        T t = (r.x() * d(0) + r.y() * d(1)) / (L * L);
        if (t < -1e-10 || t > 1 + 1e-10)
            return restore();
        t = std::min(std::max(t, T(0)), T(1));
        T hx = r.x() * nu(0) + r.y() * nu(1);
        if (std::abs(hx - crv.height(t).first) > 1e-3 * h_cell)
            return restore();
    }
    /* the curve must stay in the polygon */
    T tol = 1e-10 * h_cell;
    size_t nchk = std::max<size_t>(samples, 4 * degree);
    for (size_t k = 0; k <= nchk; k++) {
        point<T,2> x;
        vec2<T> dx;
        crv.eval(T(k) / T(nchk), x, dx);
        if (!inside(poly, halfplanes, x, tol))
            return restore();
    }

    poly.polyline.clear();
    for (size_t k = 0; k <= samples; k++) {
        point<T,2> x;
        vec2<T> dx;
        crv.eval(T(k) / T(samples), x, dx);
        poly.polyline.push_back(x);
    }
    poly.polyline.front() = poly.p0;
    poly.polyline.back() = poly.p1;

    /* each part must be star-shaped with respect to its centroid (sample
     * polygon), the center of the sectors of its quadrature */
    for (size_t side = 0; side < 2; side++) {
        auto where = (side == 0) ? location::IN_NEGATIVE_SIDE : location::IN_POSITIVE_SIDE;
        T area = 0;
        auto cg = polygon_centroid(side_polygon(poly, where), &area);
        if (area <= 0)
            return restore();
        poly.centers[side] = cg;
        if (!star_shaped(poly, crv, where, cg, std::max<size_t>(64, 8 * degree)))
            return restore();
    }

    return true;
}

/* Implicit interface: the exact zero level set, as the graph of the root of
 * the level set along the normal of the chord. Checked at "samples" points
 * (the interface must be found in the polygon, be a graph above the chord,
 * and each part must be star-shaped with respect to one of the candidate
 * centers); on success, the polyline is replaced by the samples. Returns
 * false (the polygon keeping its polyline) otherwise. */
template<typename T, typename LS2>
bool
make_implicit_interface_direct(cut_polygon<T>& poly, const LS2& phi, size_t samples)
{
    T L = (poly.p1 - poly.p0).to_vector().norm();
    if (!(L > 1e-10 * poly.h))
        return false;

    auto saved_kind = poly.kind;
    auto saved_polyline = poly.polyline;
    auto restore = [&]() {
        poly.kind = saved_kind;
        poly.polyline = saved_polyline;
        return false;
    };

    poly.kind = interface_kind::implicit;
    polygon_curve<T, LS2> crv(poly, phi);
    auto halfplanes = inward_halfplanes(poly);

    std::vector<point<T,2>> pts;
    for (size_t k = 0; k <= samples; k++) {
        point<T,2> x;
        vec2<T> dx;
        if (!crv.eval(T(k) / T(samples), x, dx))
            return restore();
        if (!inside(poly, halfplanes, x, 1e-10 * poly.h))
            return restore();
        pts.push_back(x);
    }
    pts.front() = poly.p0;
    pts.back() = poly.p1;
    poly.polyline = pts;

    for (size_t side = 0; side < 2; side++) {
        auto where = (side == 0) ? location::IN_NEGATIVE_SIDE : location::IN_POSITIVE_SIDE;
        auto tp = side_polygon(poly, where);
        T area = 0;
        std::vector<point<T,2>> candidates = {polygon_centroid(tp, &area)};
        if (area <= 0)
            return restore();
        auto chord_poly = side_polygon(poly, where, true);
        candidates.push_back(polygon_centroid(chord_poly));
        point<T,2> m;
        for (auto& p : chord_poly)
            m = m + p;
        candidates.push_back(m / T(chord_poly.size()));
        for (size_t i = 0; i < poly.vertices.size(); i++)
            if (poly.vertex_loc[i] == where)
                candidates.push_back(poly.vertices[i]);
        bool ok = false;
        for (auto& c : candidates)
            if (star_shaped(poly, crv, where, c, 64)) {
                poly.centers[side] = c;
                ok = true;
                break;
            }
        if (!ok)
            return restore();
    }

    return true;
}

template<typename T, typename LS2>
point<T,2>
find_zero_crossing(const point<T,2>& p0, const point<T,2>& p1, const LS2& phi, const T& threshold);

/* Part of a convex polygon on one side of the line m.x = d (m.x <= d if
 * keep_le): its vertices, the parent vertex of each vertex (-1 if new), and
 * the parent edge containing the edge starting at each vertex (-1 for the
 * edge on the line) */
template<typename T>
struct clipped_polygon {
    std::vector<point<T,2>> vertices;
    std::vector<int>        vertex_origin;
    std::vector<int>        edge_origin;
};

template<typename T>
clipped_polygon<T>
clip_polygon(const std::vector<point<T,2>>& v, const vec2<T>& m, T d, bool keep_le)
{
    clipped_polygon<T> raw;
    size_t n = v.size();
    std::vector<T> f(n);
    for (size_t i = 0; i < n; i++)
        f[i] = m.dot(to_vec(v[i])) - d;
    auto in = [&](size_t i) { return keep_le ? f[i] <= 0 : f[i] >= 0; };
    for (size_t i = 0; i < n; i++) {
        size_t j = (i+1) % n;
        if (in(i)) {
            raw.vertices.push_back(v[i]);
            raw.vertex_origin.push_back(int(i));
            raw.edge_origin.push_back(int(i));
        }
        if (in(i) != in(j)) {
            auto x = v[i] + (v[j] - v[i]) * (f[i] / (f[i] - f[j]));
            raw.vertices.push_back(x);
            raw.vertex_origin.push_back(-1);
            raw.edge_origin.push_back(in(i) ? -1 : int(i));
        }
    }
    /* coincident consecutive points: the last one gives the outgoing edge */
    clipped_polygon<T> ret;
    for (size_t i = 0; i < raw.vertices.size(); i++) {
        if (!ret.vertices.empty() && same_point(ret.vertices.back(), raw.vertices[i])) {
            ret.edge_origin.back() = raw.edge_origin[i];
            if (raw.vertex_origin[i] >= 0)
                ret.vertex_origin.back() = raw.vertex_origin[i];
            continue;
        }
        ret.vertices.push_back(raw.vertices[i]);
        ret.vertex_origin.push_back(raw.vertex_origin[i]);
        ret.edge_origin.push_back(raw.edge_origin[i]);
    }
    while (ret.vertices.size() > 1 && same_point(ret.vertices.front(), ret.vertices.back())) {
        if (ret.vertex_origin.back() >= 0)
            ret.vertex_origin.front() = ret.vertex_origin.back();
        ret.vertices.pop_back();
        ret.vertex_origin.pop_back();
        ret.edge_origin.pop_back();
    }
    return ret;
}

template<typename T, typename LS2>
bool
make_implicit_interface(cut_polygon<T>& poly, const LS2& phi, size_t samples = 32, size_t depth = 0,
                        size_t max_depth = 4);

/* Implicit interface of a polygon split by the line through the middle of
 * its interface, perpendicular to its chord: the parts lying on one side of
 * the interface and the cut parts, each built recursively */
template<typename T, typename LS2>
bool
split_cut_polygon(cut_polygon<T>& poly, const LS2& phi, size_t samples, size_t depth, size_t max_depth)
{
    if (!poly.convex)
        return false;

    auto saved_kind = poly.kind;
    poly.kind = interface_kind::implicit;
    polygon_curve<T, LS2> crv(poly, phi);
    point<T,2> xm = crv.chord(0.5);
    vec2<T> dxm;
    crv.eval(0.5, xm, dxm);
    poly.kind = saved_kind;

    vec2<T> tau = (poly.p1 - poly.p0).to_vector();
    tau /= tau.norm();
    T d = tau.dot(to_vec(xm));

    std::optional<point<T,2>> split_point; /* intersection of the interface with the line */
    std::vector<cut_polygon<T>> cut_parts;
    std::vector<std::pair<std::vector<point<T,2>>, location>> uncut_parts;
    for (bool keep_le : {true, false}) {
        auto cp = clip_polygon(poly.vertices, tau, d, keep_le);
        if (cp.vertices.size() < 3)
            continue;
        size_t n = cp.vertices.size();
        std::vector<location> loc(n);
        for (size_t i = 0; i < n; i++)
            loc[i] = (cp.vertex_origin[i] >= 0) ? poly.vertex_loc[cp.vertex_origin[i]]
                   : (phi(cp.vertices[i]) < 0 ? location::IN_NEGATIVE_SIDE : location::IN_POSITIVE_SIDE);
        size_t changes = 0;
        for (size_t i = 0; i < n; i++)
            if (loc[i] != loc[(i+1) % n])
                changes++;
        if (changes == 0) {
            uncut_parts.push_back(std::make_pair(cp.vertices, loc[0]));
            continue;
        }
        if (changes != 2)
            return false;
        bool ok = true;
        auto intersection = [&](size_t i) {
            int e = cp.edge_origin[i];
            if (e >= 0) {
                if (size_t(e) == poly.edge0)
                    return poly.p0;
                if (size_t(e) == poly.edge1)
                    return poly.p1;
                ok = false; /* the interface crosses an edge twice */
                return poly.p0;
            }
            if (!split_point) {
                /* the same point for both parts: from the ends ordered along tau^perp */
                auto a = cp.vertices[i], b = cp.vertices[(i+1) % n];
                vec2<T> np(-tau(1), tau(0));
                if (np.dot(to_vec(a)) > np.dot(to_vec(b)))
                    std::swap(a, b);
                split_point = find_zero_crossing(a, b, phi, poly.h / 1e20);
            }
            return *split_point;
        };
        cut_polygon<T> part;
        if (!init_cut_polygon(part, cp.vertices, loc, std::vector<size_t>(), intersection) || !ok)
            return false;
        if (!make_implicit_interface(part, phi, samples, depth + 1, max_depth))
            return false;
        cut_parts.push_back(part);
    }
    if (cut_parts.empty())
        return false;

    /* order the cut parts along the interface */
    std::vector<cut_polygon<T>> ordered;
    point<T,2> current = poly.p0;
    std::vector<bool> used(cut_parts.size(), false);
    for (size_t k = 0; k < cut_parts.size(); k++) {
        bool found = false;
        for (size_t j = 0; j < cut_parts.size() && !found; j++) {
            if (used[j] || !same_point(cut_parts[j].p0, current))
                continue;
            ordered.push_back(cut_parts[j]);
            current = cut_parts[j].p1;
            used[j] = found = true;
        }
        if (!found)
            return false;
    }
    if (!same_point(current, poly.p1))
        return false;

    poly.children = ordered;
    poly.uncut_children = uncut_parts;
    poly.kind = interface_kind::implicit;
    poly.polyline.clear();
    for (auto& c : poly.children)
        for (auto& p : c.polyline)
            if (poly.polyline.empty() || !same_point(poly.polyline.back(), p))
                poly.polyline.push_back(p);
    return true;
}

/* Implicit interface (see make_implicit_interface_direct). If a part of the
 * polygon is not star-shaped, the polygon is split (up to max_depth times).
 * Returns false (the polygon keeping its polyline) on failure. */
template<typename T, typename LS2>
bool
make_implicit_interface(cut_polygon<T>& poly, const LS2& phi, size_t samples, size_t depth, size_t max_depth)
{
    poly.children.clear();
    poly.uncut_children.clear();
    if (make_implicit_interface_direct(poly, phi, samples))
        return true;
    if (depth >= max_depth)
        return false;
    auto saved_kind = poly.kind;
    auto saved_polyline = poly.polyline;
    if (split_cut_polygon(poly, phi, samples, depth, max_depth))
        return true;
    poly.children.clear();
    poly.uncut_children.clear();
    poly.kind = saved_kind;
    poly.polyline = saved_polyline;
    return false;
}

/* The leaves of the decomposition of an implicit interface, ordered along
 * the interface */
template<typename T>
void
interface_leaves(const cut_polygon<T>& poly, std::vector<const cut_polygon<T>*>& out)
{
    if (poly.children.empty()) {
        out.push_back(&poly);
        return;
    }
    for (auto& c : poly.children)
        interface_leaves(c, out);
}

/* Zero of the level set on the segment [p0, p1] (bisection, ProtoN) */
template<typename T, typename LS2>
point<T,2>
find_zero_crossing(const point<T,2>& p0, const point<T,2>& p1, const LS2& phi, const T& threshold) {

    /* !!! We assume that the level set function *has* a zero crossing
     * between p0 and p1 !!! */
    auto pa = p0;
    auto pb = p1;
    auto pm = (pa+pb)/2.0;
    auto pm_prev = pm;

    T x_diff_sq, y_diff_sq;
    size_t max_iter = 50;

    do {
        auto lb = phi(pb);
        auto lm = phi(pm);

        if ( (lb >= 0 && lm >= 0) || (lb < 0 && lm < 0) ) {   /* intersection is between pa and pm */
            pm_prev = pm;
            pb = pm;
            pm = (pa+pb)/2.0;
        }
        else {   /* intersection is between pm and pb */
            pm_prev = pm;
            pa = pm;
            pm = (pa+pb)/2.0;
        }

        x_diff_sq = (pm_prev.x() - pm.x()) * (pm_prev.x() - pm.x());
        y_diff_sq = (pm_prev.y() - pm.y()) * (pm_prev.y() - pm.y());

    } while ( (std::sqrt(x_diff_sq + y_diff_sq) > threshold) && max_iter-- );

    return pm;
}

/* Recursive splitting of the interface between the points min and max */
template<typename T, typename LS2>
void
refine_polyline(cut_polygon<T>& poly, const LS2& phi, size_t min, size_t max)
{
    if ((max-min) < 2)
        return;

    size_t mid = (max+min)/2;
    auto& itf = poly.polyline;
    auto p0 = itf.at(min);
    auto p1 = itf.at(max);
    auto pm = (p0+p1)/2.0;
    auto pt = p1 - p0;
    auto pn = point<T,2>(-pt.y(), pt.x());
    auto ps1 = pm + pn;
    auto ps2 = pm - pn;

    auto lm = phi(pm);
    auto ls1 = phi(ps1);
    auto ls2 = phi(ps2);

    point<T,2> ip;
    auto threshold = poly.h / 1e20;
    if ( !((lm >= 0 && ls1 >= 0) || (lm < 0 && ls1 < 0)) )
        ip = find_zero_crossing(pm, ps1, phi, threshold);
    else if ( !((lm >= 0 && ls2 >= 0) || (lm < 0 && ls2 < 0)) )
        ip = find_zero_crossing(pm, ps2, phi, threshold);
    else
        throw std::logic_error("interface not found in search range: the mesh does not resolve the interface");

    itf.at(mid) = ip;

    refine_polyline(poly, phi, min, mid);
    refine_polyline(poly, phi, mid, max);
}

/* Polyline of 2^levels segments (ProtoN's refine_interface). The points
 * outside the polygon are projected onto its boundary: the interface then
 * crosses an edge twice, and the small part of the subdomain lying in the
 * neighboring element is neglected (an element is cut once). Returns true if
 * a point was projected. */
template<typename T, typename LS2>
bool
refine_polyline(cut_polygon<T>& poly, const LS2& phi, size_t levels)
{
    poly.kind = interface_kind::polyline;
    if (levels == 0) {
        poly.polyline = {poly.p0, poly.p1};
        return false;
    }

    size_t npts = size_t(1) << levels;
    poly.polyline.assign(npts + 1, point<T,2>());
    poly.polyline.front() = poly.p0;
    poly.polyline.back()  = poly.p1;
    refine_polyline(poly, phi, 0, npts);

    /* axis-aligned rectangles: clamping to the box, as in ProtoN */
    T b[4] = {poly.vertices[0].x(), poly.vertices[0].x(), poly.vertices[0].y(), poly.vertices[0].y()};
    for (auto& v : poly.vertices) {
        b[0] = std::min(b[0], v.x()); b[1] = std::max(b[1], v.x());
        b[2] = std::min(b[2], v.y()); b[3] = std::max(b[3], v.y());
    }
    T tol_r = 1e-12 * std::max(b[1] - b[0], b[3] - b[2]);
    bool rectangle = (poly.vertices.size() == 4);
    for (auto& v : poly.vertices) {
        bool cx = std::abs(v.x() - b[0]) <= tol_r || std::abs(v.x() - b[1]) <= tol_r;
        bool cy = std::abs(v.y() - b[2]) <= tol_r || std::abs(v.y() - b[3]) <= tol_r;
        rectangle = rectangle && cx && cy;
    }

    T tol = rectangle ? 1e-12 * (b[1] - b[0]) : 1e-12 * poly.h;
    auto halfplanes = inward_halfplanes(poly);
    bool clamped = false;
    for (auto& p : poly.polyline) {
        if (rectangle) {
            auto x = std::min(std::max(p.x(), b[0]), b[1]);
            auto y = std::min(std::max(p.y(), b[2]), b[3]);
            if (std::abs(x - p.x()) > tol || std::abs(y - p.y()) > tol) {
                p = point<T,2>(x, y);
                clamped = true;
            }
            continue;
        }
        if (!poly.convex || inside(poly, halfplanes, p, tol))
            continue;
        /* closest point of the boundary */
        point<T,2> best = p;
        T best_d = std::numeric_limits<T>::max();
        size_t n = poly.vertices.size();
        for (size_t i = 0; i < n; i++) {
            const auto& a = poly.vertices[i];
            const auto& bb = poly.vertices[(i+1)%n];
            auto e = (bb - a).to_vector();
            T l2 = e.squaredNorm();
            T t = (l2 > 0) ? std::min(std::max((p - a).to_vector().dot(e) / l2, T(0)), T(1)) : T(0);
            auto q = a + (bb - a) * t;
            T dq = (p - q).to_vector().norm();
            if (dq < best_d) {
                best_d = dq;
                best = q;
            }
        }
        p = best;
        clamped = true;
    }
    if (clamped) {
        /* remove the points made coincident by the projection */
        std::vector<point<T,2>> pts;
        for (auto& p : poly.polyline)
            if (pts.empty() || (p - pts.back()).to_vector().norm() > tol)
                pts.push_back(p);
        poly.polyline = pts;
    }

    return clamped;
}

/* Triangulation of a simple counterclockwise polygon by ear clipping */
template<typename T>
std::vector<std::array<point<T,2>, 3>>
ear_clipping(const std::vector<point<T,2>>& polygon)
{
    std::vector<std::array<point<T,2>, 3>> ret;
    std::vector<point<T,2>> pts;
    for (auto& p : polygon)
        if (pts.empty() || !same_point(p, pts.back()))
            pts.push_back(p);
    while (pts.size() > 1 && same_point(pts.front(), pts.back()))
        pts.pop_back();

    T scale = 0;
    for (auto& p : pts)
        scale = std::max(scale, (p - pts[0]).to_vector().norm());
    T tol = 1e-14 * scale * scale;

    auto in_triangle = [&](const point<T,2>& p, const point<T,2>& a, const point<T,2>& b, const point<T,2>& c) {
        return cross(b - a, p - a) > tol && cross(c - b, p - b) > tol && cross(a - c, p - c) > tol;
    };

    while (pts.size() > 3) {
        size_t n = pts.size();
        bool clipped = false;
        for (size_t i = 0; i < n && !clipped; i++) {
            const auto& a = pts[(i + n - 1) % n];
            const auto& b = pts[i];
            const auto& c = pts[(i + 1) % n];
            if (cross(b - a, c - b) < -tol)
                continue; /* reflex vertex */
            bool ear = true;
            for (size_t j = 0; j < n && ear; j++) {
                if (j == i || j == (i + n - 1) % n || j == (i + 1) % n)
                    continue;
                if (in_triangle(pts[j], a, b, c))
                    ear = false;
            }
            if (!ear)
                continue;
            ret.push_back({a, b, c});
            pts.erase(pts.begin() + i);
            clipped = true;
        }
        if (!clipped) {
            /* no ear (round-off): fan of the remaining polygon */
            for (size_t i = 1; i + 1 < pts.size(); i++)
                ret.push_back({pts[0], pts[i], pts[i+1]});
            return ret;
        }
    }
    if (pts.size() == 3)
        ret.push_back({pts[0], pts[1], pts[2]});
    return ret;
}

/* Quadrature of degree "degree" on the triangle (a, b, c) (signed weights
 * if the triangle is clockwise) */
template<typename T>
void
add_triangle_quadrature(std::vector<std::pair<point<T,2>, T>>& ret, size_t degree,
                        const point<T,2>& a, const point<T,2>& b, const point<T,2>& c)
{
    auto v1 = b - a;
    auto v2 = c - a;
    T det = cross(v1, v2);
    for (auto& qp : reference_triangle_rule<T>(degree))
        ret.push_back(std::make_pair(a + v1*qp.first.x() + v2*qp.first.y(), qp.second*det));
}

/* Quadrature of degree "degree" on the part "where" of a cut polygon. With
 * the polynomial and implicit interfaces, the part is split into the
 * triangles (c, straight edge) and the curved sector X(xi, t) = c + xi
 * (x(t) - c), of Jacobian xi det(x(t) - c, x'(t)), c being its center: all
 * the points lie in the part and all the weights are positive. The polyline
 * is triangulated by ear clipping. */
template<typename T, typename LS2>
std::vector<std::pair<point<T,2>, T>>
side_quadrature(const cut_polygon<T>& poly, const LS2& phi, size_t degree, location where)
{
    std::vector<std::pair<point<T,2>, T>> ret;
    bool neg = (where == location::IN_NEGATIVE_SIDE);

    if (!poly.children.empty()) {
        for (auto& c : poly.children) {
            auto q = side_quadrature(c, phi, degree, where);
            ret.insert(ret.end(), q.begin(), q.end());
        }
        for (auto& u : poly.uncut_children) {
            if (u.second != where)
                continue;
            for (size_t i = 1; i + 1 < u.first.size(); i++)
                add_triangle_quadrature(ret, degree, u.first[0], u.first[i], u.first[i+1]);
        }
        return ret;
    }

    if (poly.kind == interface_kind::polyline) {
        for (auto& tri : ear_clipping(side_polygon(poly, where)))
            add_triangle_quadrature(ret, degree, tri[0], tri[1], tri[2]);
        return ret;
    }

    const auto& c = poly.centers[neg ? 0 : 1];
    auto chord_poly = side_polygon(poly, where, true);
    for (size_t i = 0; i < chord_poly.size(); i++) {
        const auto& a = chord_poly[i];
        const auto& b = chord_poly[(i+1) % chord_poly.size()];
        if (is_chord(poly, a, b, where))
            continue;
        add_triangle_quadrature(ret, degree, c, a, b);
    }

    polygon_curve<T, LS2> crv(poly, phi);
    const auto& rule_xi = gauss_legendre_01<T>((degree + 1) / 2 + 1);
    for (auto& sw : gauss_legendre_01<T>(crv.nodes_count(degree, true))) {
        T t = neg ? sw.first : 1.0 - sw.first;
        point<T,2> x = crv.chord(t);
        vec2<T> dx = (poly.p1 - poly.p0).to_vector();
        crv.eval(t, x, dx);
        auto g = x - c;
        T det = cross(g, dx);
        if (!neg)
            det = -det;
        for (auto& xw : rule_xi)
            ret.push_back(std::make_pair(c + g * xw.first, sw.second * xw.second * xw.first * det));
    }
    return ret;
}

/* Quadrature on the interface of a cut polygon: points, weights (length)
 * and unit normals, pointing from the negative to the positive side */
template<typename T>
struct curve_quadrature_t {
    std::vector<std::pair<point<T,2>, T>> qps;
    std::vector<vec2<T>>                   normals;
};

template<typename T, typename LS2>
curve_quadrature_t<T>
curve_quadrature(const cut_polygon<T>& poly, const LS2& phi, size_t degree)
{
    curve_quadrature_t<T> ret;
    if (!poly.children.empty()) {
        for (auto& c : poly.children) {
            auto q = curve_quadrature(c, phi, degree);
            ret.qps.insert(ret.qps.end(), q.qps.begin(), q.qps.end());
            ret.normals.insert(ret.normals.end(), q.normals.begin(), q.normals.end());
        }
        return ret;
    }
    if (poly.kind == interface_kind::polyline) {
        const auto& itf = poly.polyline;
        for (size_t i = 1; i < itf.size(); i++) {
            auto p0 = itf[i-1], p1 = itf[i];
            auto e = (p1 - p0).to_vector();
            T meas = e.norm();
            vec2<T> n(e(1)/meas, -e(0)/meas);
            for (auto& qp : disk::quadrature::gauss_legendre(degree, p0, p1)) {
                ret.qps.push_back(std::make_pair(qp.point(), qp.weight()));
                ret.normals.push_back(n);
            }
        }
        return ret;
    }

    polygon_curve<T, LS2> crv(poly, phi);
    for (auto& tw : gauss_legendre_01<T>(crv.nodes_count(degree, false))) {
        point<T,2> x = crv.chord(tw.first);
        vec2<T> dx = (poly.p1 - poly.p0).to_vector();
        crv.eval(tw.first, x, dx);
        T ds = dx.norm();
        ret.qps.push_back(std::make_pair(x, tw.second * ds));
        ret.normals.push_back(vec2<T>(dx(1), -dx(0)) / ds);
    }
    return ret;
}

/***************************************************************************/
/*                 TETRAHEDRA CUT BY THE INTERFACE (FALLBACK)               */
/***************************************************************************/

/* Quadratures of the two parts of a polytope and of the interface in it */
template<typename T, size_t DIM>
struct cut_quadratures {
    std::vector<quadrature_point<T, DIM>>   neg, pos;   // parts
    std::vector<quadrature_point<T, DIM>>   itf;        // interface
    std::vector<Eigen::Matrix<T, DIM, 1>>   normals;    // from the negative to the positive side
};

template<typename T>
void
add_tetrahedron_quadrature(std::vector<quadrature_point<T,3>>& ret, size_t degree,
                           const point<T,3>& a, const point<T,3>& b, const point<T,3>& c, const point<T,3>& d)
{
    Eigen::Matrix<T,3,3> J;
    J.col(0) = (b - a).to_vector();
    J.col(1) = (c - a).to_vector();
    J.col(2) = (d - a).to_vector();
    T det = std::abs(J.determinant());
    for (auto& qp : reference_tetrahedron_rule<T>(degree)) {
        Eigen::Matrix<T,3,1> x = a.to_vector() + J * Eigen::Matrix<T,3,1>(qp.first[0], qp.first[1], qp.first[2]);
        ret.push_back(make_qp(point<T,3>(x), qp.second * det));
    }
}

template<typename T>
void
add_triangle_quadrature(std::vector<quadrature_point<T,3>>& ret, size_t degree,
                        const point<T,3>& a, const point<T,3>& b, const point<T,3>& c)
{
    auto v1 = (b - a).to_vector();
    auto v2 = (c - a).to_vector();
    T area2 = v1.cross(v2).norm();
    for (auto& qp : reference_triangle_rule<T>(degree)) {
        Eigen::Matrix<T,3,1> x = a.to_vector() + v1*qp.first.x() + v2*qp.first.y();
        ret.push_back(make_qp(point<T,3>(x), qp.second * area2));
    }
}

/* Quadratures of the parts of the tetrahedron "tet" (degree "degree"),
 * refined "levels" times (8 subtetrahedra, Bey's splitting): in each final
 * tetrahedron, the interface is approximated by the plane polygon through
 * its intersection points with the edges (exact roots of the level set). */
template<typename T, typename LS>
void
subdivide_tetrahedron(const std::array<point<T,3>, 4>& tet, const LS& ls, size_t levels, size_t degree,
                      cut_quadratures<T,3>& out)
{
    if (levels > 0) {
        std::array<T,4> f;
        for (size_t i = 0; i < 4; i++)
            f[i] = ls(tet[i]);
        bool all_neg = std::all_of(f.begin(), f.end(), [](T v) { return v < 0; });
        bool all_pos = std::all_of(f.begin(), f.end(), [](T v) { return v >= 0; });
        /* the interface may still cross a tetrahedron whose vertices are on
         * one side; the subdivision is stopped only if the level set does not
         * vanish in a ball containing the tetrahedron, roughly */
        point<T,3> g = (tet[0] + tet[1] + tet[2] + tet[3]) / T(4);
        T r = 0;
        for (auto& p : tet)
            r = std::max(r, (p - g).to_vector().norm());
        T gn = ls.gradient(g).norm();
        bool far = std::abs(ls(g)) > 2.0 * gn * r;
        if ((all_neg || all_pos) && far) {
            add_tetrahedron_quadrature(all_neg ? out.neg : out.pos, degree, tet[0], tet[1], tet[2], tet[3]);
            return;
        }
        auto m = [&](size_t i, size_t j) { return (tet[i] + tet[j]) / T(2); };
        auto p01 = m(0,1), p02 = m(0,2), p03 = m(0,3), p12 = m(1,2), p13 = m(1,3), p23 = m(2,3);
        std::array<std::array<point<T,3>,4>, 8> sub = {{
            {tet[0], p01, p02, p03}, {p01, tet[1], p12, p13}, {p02, p12, tet[2], p23}, {p03, p13, p23, tet[3]},
            {p01, p02, p03, p13}, {p01, p02, p12, p13}, {p02, p03, p13, p23}, {p02, p12, p13, p23}
        }};
        for (auto& s : sub)
            subdivide_tetrahedron(s, ls, levels - 1, degree, out);
        return;
    }

    std::array<T,4> f;
    std::array<bool,4> neg;
    size_t nneg = 0;
    for (size_t i = 0; i < 4; i++) {
        f[i] = ls(tet[i]);
        neg[i] = (f[i] < 0);
        nneg += neg[i] ? 1 : 0;
    }
    if (nneg == 0 || nneg == 4) {
        add_tetrahedron_quadrature(nneg == 4 ? out.neg : out.pos, degree, tet[0], tet[1], tet[2], tet[3]);
        return;
    }

    /* root of the level set on the edge (i, j) */
    auto q = [&](size_t i, size_t j) {
        auto a = tet[i], b = tet[j];
        auto e = b - a;
        auto fdf = [&](T s) {
            auto x = a + e * s;
            return std::make_pair(ls(x), ls.gradient(x).dot(e.to_vector()));
        };
        auto s = find_root_near(fdf, T(0), T(1), T(0.5), T(1e-3), T(1e-15));
        return a + e * (s ? *s : f[i]/(f[i] - f[j]));
    };

    /* orient an interface triangle from the negative to the positive side */
    auto add_interface = [&](const point<T,3>& a, const point<T,3>& b, const point<T,3>& c, const point<T,3>& xneg) {
        Eigen::Matrix<T,3,1> n = (b - a).to_vector().cross((c - a).to_vector());
        T nn = n.norm();
        if (nn == 0)
            return;
        n /= nn;
        if (n.dot((xneg - a).to_vector()) > 0)
            n = -n;
        size_t before = out.itf.size();
        add_triangle_quadrature(out.itf, degree, a, b, c);
        for (size_t k = before; k < out.itf.size(); k++)
            out.normals.push_back(n);
    };

    /* prism (a, b, c | a', b', c'), a' facing a */
    auto add_prism = [&](std::vector<quadrature_point<T,3>>& dst, const point<T,3>& a, const point<T,3>& b,
                         const point<T,3>& c, const point<T,3>& a2, const point<T,3>& b2, const point<T,3>& c2) {
        add_tetrahedron_quadrature(dst, degree, a, b, c, a2);
        add_tetrahedron_quadrature(dst, degree, b, c, a2, b2);
        add_tetrahedron_quadrature(dst, degree, c, a2, b2, c2);
    };

    if (nneg == 1 || nneg == 3) {
        /* vertex i alone on its side */
        size_t i = 0;
        for (size_t k = 0; k < 4; k++)
            if ((nneg == 1) == neg[k])
                i = k;
        std::array<size_t,3> o;
        size_t c = 0;
        for (size_t k = 0; k < 4; k++)
            if (k != i)
                o[c++] = k;
        auto qa = q(i, o[0]), qb = q(i, o[1]), qc = q(i, o[2]);
        auto& alone = neg[i] ? out.neg : out.pos;
        auto& others = neg[i] ? out.pos : out.neg;
        add_tetrahedron_quadrature(alone, degree, tet[i], qa, qb, qc);
        add_prism(others, tet[o[0]], tet[o[1]], tet[o[2]], qa, qb, qc);
        add_interface(qa, qb, qc, neg[i] ? tet[i] : tet[o[0]]);
        return;
    }

    /* two vertices on each side: (i, j) negative, (k, l) positive */
    std::array<size_t,2> a, b;
    size_t na = 0, nb = 0;
    for (size_t k = 0; k < 4; k++) {
        if (neg[k])
            a[na++] = k;
        else
            b[nb++] = k;
    }
    auto q_ik = q(a[0], b[0]), q_il = q(a[0], b[1]), q_jk = q(a[1], b[0]), q_jl = q(a[1], b[1]);
    add_prism(out.neg, tet[a[0]], q_ik, q_il, tet[a[1]], q_jk, q_jl);
    add_prism(out.pos, tet[b[0]], q_ik, q_jk, tet[b[1]], q_il, q_jl);
    add_interface(q_ik, q_il, q_jl, tet[a[0]]);
    add_interface(q_ik, q_jl, q_jk, tet[a[0]]);
}

/***************************************************************************/
/*                      POLYHEDRA CUT BY THE INTERFACE                      */
/***************************************************************************/

/* Unit normal of a planar polygon (Newell's method): the vertices are
 * counterclockwise around it */
template<typename Points>
Eigen::Matrix<typename Points::value_type::value_type, 3, 1>
newell_normal(const Points& pts)
{
    using T = typename Points::value_type::value_type;
    Eigen::Matrix<T,3,1> n = Eigen::Matrix<T,3,1>::Zero();
    for (size_t i = 0; i < pts.size(); i++) {
        const auto& a = pts[i];
        const auto& b = pts[(i+1) % pts.size()];
        n(0) += (a.y() - b.y()) * (a.z() + b.z());
        n(1) += (a.z() - b.z()) * (a.x() + b.x());
        n(2) += (a.x() - b.x()) * (a.y() + b.y());
    }
    return n / n.norm();
}

/* Orthonormal frame (origin, U) of the plane of a polygon of space, in which
 * its vertices are counterclockwise */
template<typename Points, typename T>
void
polygon_frame(const Points& pts, point<T,3>& origin, Eigen::Matrix<T,3,2>& U)
{
    origin = point<T,3>();
    for (auto& p : pts)
        origin = origin + p;
    origin = origin / T(pts.size());
    auto n = newell_normal(pts);
    Eigen::Matrix<T,3,1> e1 = (pts[1] - pts[0]).to_vector();
    e1 -= e1.dot(n) * n;
    e1 /= e1.norm();
    Eigen::Matrix<T,3,1> e2 = n.cross(e1);
    U.col(0) = e1;
    U.col(1) = e2;
}

/* lexicographic order of two points */
template<typename T, size_t DIM>
bool
point_less(const point<T,DIM>& a, const point<T,DIM>& b)
{
    for (size_t i = 0; i < DIM; i++) {
        if (a[i] < b[i])
            return true;
        if (a[i] > b[i])
            return false;
    }
    return false;
}

/* An edge of a polyhedron, by the coordinates of its ends */
template<typename T>
using edge_key = std::array<T,6>;

template<typename T>
edge_key<T>
make_edge_key(const point<T,3>& a, const point<T,3>& b)
{
    if (point_less(b, a))
        return {b.x(), b.y(), b.z(), a.x(), a.y(), a.z()};
    return {a.x(), a.y(), a.z(), b.x(), b.y(), b.z()};
}

/* Zero of the level set on the edge [a, b] (bisection from the smaller end:
 * both faces of the edge get the same point) */
template<typename T, typename LS>
point<T,3>
edge_zero(point<T,3> a, point<T,3> b, const LS& ls)
{
    if (point_less(b, a))
        std::swap(a, b);
    auto fa = ls(a);
    for (size_t it = 0; it < 60; it++) {
        auto pm = (a + b) / 2.0;
        auto fm = ls(pm);
        if ((fa < 0) == (fm < 0)) {
            a = pm;
            fa = fm;
        }
        else
            b = pm;
    }
    return (a + b) / 2.0;
}

/* Quadrature of degree "degree" on a convex polygon of space (positive
 * weights): tensor Gauss rule on a parallelogram, fan of triangles from a
 * vertex otherwise */
template<typename T>
void
add_polygon_quadrature(std::vector<quadrature_point<T,3>>& ret, size_t degree, const std::vector<point<T,3>>& v)
{
    if (v.size() == 4) {
        vec3<T> e1 = (v[1] - v[0]).to_vector(), e3 = (v[3] - v[0]).to_vector();
        vec3<T> diag = (v[2] - v[0]).to_vector();
        if ((diag - e1 - e3).norm() <= 1e-12 * diag.norm()) {
            T area = e1.cross(e3).norm();
            const auto& g = gauss_legendre_01<T>(degree/2 + 1);
            for (auto& gs : g)
                for (auto& gt : g)
                    ret.push_back(make_qp(v[0] + point<T,3>(vec3<T>(e1 * gs.first + e3 * gt.first)),
                                          gs.second * gt.second * area));
            return;
        }
    }
    for (size_t i = 1; i + 1 < v.size(); i++)
        add_triangle_quadrature(ret, degree, v[0], v[i], v[i+1]);
}

/* Interface curve of a cut face of a polyhedron, in the coordinates of space:
 * the leaves of the decomposition of its implicit interface, or the segments
 * of its polyline, from p0 to p1 */
template<typename T>
class face_interface_curve {

    const cut_polyhedron_face<T>&       m_f;
    planar_level_set<T,3>               m_phi;
    std::vector<const cut_polygon<T>*>  m_leaves;

public:
    face_interface_curve(const level_set<T,3>& ls, const cut_polyhedron_face<T>& f)
        : m_f(f), m_phi(ls, f.origin, f.U)
    {
        if (curved())
            interface_leaves(f.poly, m_leaves);
    }

    face_interface_curve(const face_interface_curve&) = delete;

    bool curved() const { return m_f.poly.kind == interface_kind::implicit; }

    size_t pieces() const { return curved() ? m_leaves.size() : m_f.poly.polyline.size() - 1; }

    /* point and derivative with respect to t in [0,1] on the piece k */
    void
    eval(size_t k, T t, point<T,3>& x, Eigen::Matrix<T,3,1>& dx) const
    {
        if (!curved()) {
            const auto& a = m_f.poly.polyline[k];
            const auto& b = m_f.poly.polyline[k+1];
            x  = m_phi.to_global(a + (b - a) * t);
            dx = m_phi.to_global_vector((b - a).to_vector());
            return;
        }
        const auto& leaf = *m_leaves[k];
        polygon_curve<T, planar_level_set<T,3>> crv(leaf, m_phi);
        point<T,2> u = crv.chord(t);
        vec2<T> du = (leaf.p1 - leaf.p0).to_vector();
        crv.eval(t, u, du);
        x  = m_phi.to_global(u);
        dx = m_phi.to_global_vector(du);
    }
};

/* Cut of a face of a polyhedron (vertices counterclockwise seen from
 * outside): frame of its plane, cut polygon with the implicit interface (the
 * polyline if it fails). Returns false if the face is not cut exactly once. */
template<typename T, typename LS>
bool
make_cut_face(cut_polyhedron_face<T>& f, const std::vector<location>& loc, const LS& ls,
              std::map<edge_key<T>, point<T,3>>& edge_points, size_t refine_levels)
{
    polygon_frame(f.vertices, f.origin, f.U);
    planar_level_set<T,3> phi(ls, f.origin, f.U);
    size_t n = f.vertices.size();
    std::vector<point<T,2>> lv;
    for (auto& v : f.vertices)
        lv.push_back(phi.to_local(v));
    auto intersection = [&](size_t i) {
        const auto& a = f.vertices[i];
        const auto& b = f.vertices[(i+1) % n];
        auto key = make_edge_key(a, b);
        auto it = edge_points.find(key);
        if (it == edge_points.end())
            it = edge_points.insert(std::make_pair(key, edge_zero(a, b, ls))).first;
        return phi.to_local(it->second);
    };
    if (!init_cut_polygon(f.poly, lv, loc, std::vector<size_t>(), intersection))
        return false;
    if (!make_implicit_interface(f.poly, phi))
        refine_polyline(f.poly, phi, std::max<size_t>(refine_levels, 1));
    return true;
}

/* Lift of the point y of the plane of a cut polyhedron onto the interface:
 * root of the level set along nu, closest to the height s0 */
template<typename T, typename LS>
std::optional<point<T,3>>
lift_to_interface(const LS& ls, const point<T,3>& y, const vec3<T>& nu, T s0, T h)
{
    auto fdf = [&](T s) {
        point<T,3> x = y + point<T,3>(vec3<T>(nu * s));
        return std::make_pair(ls(x), ls.gradient(x).dot(nu));
    };
    auto s = find_root_near(fdf, -h, h, s0, T(1e-3) * h, T(1e-15) * h);
    if (!s)
        return std::nullopt;
    return y + point<T,3>(vec3<T>(nu * (*s)));
}

/* Interface patch of a cut polyhedron whose faces are cut: loop of their
 * interface curves, plane and region D, centers of the cones of the two
 * parts. Returns false if the polyhedron must be split. */
template<typename T, typename LS>
bool
make_interface_patch(cut_polyhedron<T>& P, const LS& ls)
{
    using point3 = point<T,3>;
    T h = P.h;
    T tol = 1e-10 * h;

    /* cut faces and their cut edges */
    struct cut_face { size_t idx; edge_key<T> e0, e1; };
    std::vector<cut_face> cfs;
    for (size_t i = 0; i < P.faces.size(); i++) {
        const auto& f = P.faces[i];
        if (f.loc != location::ON_INTERFACE)
            continue;
        size_t n = f.vertices.size();
        auto key = [&](size_t e) { return make_edge_key(f.vertices[e], f.vertices[(e+1) % n]); };
        cfs.push_back({i, key(f.poly.edge0), key(f.poly.edge1)});
    }
    if (cfs.size() < 3)
        return false;

    /* loop: each cut edge is shared by two cut faces */
    std::vector<bool> used(cfs.size(), false);
    P.loop.clear();
    P.loop.push_back(std::make_pair(cfs[0].idx, false));
    used[0] = true;
    auto exit_edge = cfs[0].e1;
    auto start_edge = cfs[0].e0;
    while (exit_edge != start_edge) {
        bool found = false;
        for (size_t j = 0; j < cfs.size() && !found; j++) {
            if (used[j])
                continue;
            if (cfs[j].e0 == exit_edge) {
                P.loop.push_back(std::make_pair(cfs[j].idx, false));
                exit_edge = cfs[j].e1;
                used[j] = found = true;
            }
            else if (cfs[j].e1 == exit_edge) {
                P.loop.push_back(std::make_pair(cfs[j].idx, true));
                exit_edge = cfs[j].e0;
                used[j] = found = true;
            }
        }
        if (!found)
            return false;
    }
    if (P.loop.size() != cfs.size()) {
        P.cut_more_than_once = true;
        return false;
    }

    /* loop points (entries of the faces), plane, orientation */
    std::vector<point3> lp;
    for (auto& e : P.loop) {
        const auto& f = P.faces[e.first];
        planar_level_set<T,3> phi(ls, f.origin, f.U);
        lp.push_back(phi.to_global(e.second ? f.poly.p1 : f.poly.p0));
    }
    point3 o;
    vec3<T> nu = vec3<T>::Zero();
    for (auto& p : lp) {
        o = o + p;
        vec3<T> g = ls.gradient(p);
        nu += g / g.norm();
    }
    o = o / T(lp.size());
    if (nu.norm() < 1e-3 * T(lp.size()))
        return false;
    nu /= nu.norm();
    size_t ax = 0;
    for (size_t a = 1; a < 3; a++)
        if (std::abs(nu(a)) < std::abs(nu(ax)))
            ax = a;
    vec3<T> ea = vec3<T>::Zero();
    ea(ax) = 1;
    vec3<T> a1 = nu.cross(ea);
    a1 /= a1.norm();
    vec3<T> a2 = nu.cross(a1);
    P.plane_origin = o;
    P.nu = nu;
    P.a1 = a1;
    P.a2 = a2;

    auto project = [&](const point3& x) {
        vec3<T> r = (x - o).to_vector();
        return point<T,2>(r.dot(a1), r.dot(a2));
    };
    std::vector<point<T,2>> lp2;
    for (auto& p : lp)
        lp2.push_back(project(p));
    T area = 0;
    auto cD = polygon_centroid(lp2, &area);
    if (std::abs(area) < 1e-12 * h * h)
        return false;
    if (area < 0) {
        /* counterclockwise around nu */
        std::reverse(P.loop.begin(), P.loop.end());
        for (auto& e : P.loop)
            e.second = !e.second;
    }
    P.D_center = cD;
    auto yc = o + point3(vec3<T>(a1 * cD.x() + a2 * cD.y()));
    auto xc = lift_to_interface(ls, yc, nu, T(0), h);
    if (!xc)
        return false;
    P.center_height = ((*xc) - yc).to_vector().dot(nu);

    /* D star-shaped with respect to its center, interface graph above D:
     * checked on samples of the boundary and of the sectors */
    P.interface.clear();
    P.interface.push_back(*xc);
    for (auto& entry : P.loop) {
        face_interface_curve<T> fcrv(ls, P.faces[entry.first]);
        size_t np = fcrv.pieces();
        size_t ns = fcrv.curved() ? 16 : 2;
        for (size_t k = 0; k < np; k++) {
            for (size_t j = 0; j <= ns; j++) {
                T t = T(j) / T(ns);
                point3 X;
                vec3<T> dX;
                fcrv.eval(k, t, X, dX);
                if (entry.second)
                    dX = -dX;
                auto b = project(X);
                vec2<T> db(dX.dot(a1), dX.dot(a2));
                if (cross(b - cD, db) < -1e-12 * h * h)
                    return false;
                T sb = (X - o).to_vector().dot(nu);
                for (T xi : {T(0.25), T(0.5), T(0.75)}) {
                    auto y = cD + (b - cD) * xi;
                    auto x = lift_to_interface(ls, o + point3(vec3<T>(a1 * y.x() + a2 * y.y())), nu,
                                               (1 - xi) * P.center_height + xi * sb, h);
                    if (!x)
                        return false;
                    vec3<T> g = ls.gradient(*x);
                    if (g.dot(nu) <= 0.2 * g.norm())
                        return false;
                    if (xi == T(0.5) && j % 4 == 0)
                        P.interface.push_back(*x);
                }
                if (j % 4 == 0)
                    P.interface.push_back(X);
            }
        }
    }

    /* the interface must lie in the polyhedron (half-spaces of its faces) */
    for (auto& x : P.interface)
        for (auto& f : P.faces)
            if (f.normal.dot((x - f.vertices[0]).to_vector()) > tol)
                return false;

    /* centers of the cones of the two parts */
    std::vector<point3> vertices[2];
    for (auto& f : P.faces)
        for (auto& v : f.vertices) {
            size_t side = (ls(v) < 0) ? 0 : 1;
            bool known = false;
            for (auto& w : vertices[side])
                known = known || same_point(v, w);
            if (!known)
                vertices[side].push_back(v);
        }
    for (size_t side = 0; side < 2; side++) {
        auto where = (side == 0) ? location::IN_NEGATIVE_SIDE : location::IN_POSITIVE_SIDE;
        T sgn = (side == 0) ? 1.0 : -1.0; /* outward normal of the part on the interface: sgn grad phi */
        const auto& vs = vertices[side];
        if (vs.empty())
            return false;
        std::vector<point3> candidates;
        point3 m;
        for (auto& v : vs)
            m = m + v;
        point3 mv = m / T(vs.size());
        for (auto& p : lp)
            m = m + p;
        candidates.push_back(m / T(vs.size() + lp.size()));
        candidates.push_back(mv);
        for (T a : {T(0.1), T(0.25), T(0.4)})
            candidates.push_back(o - point3(vec3<T>(nu * (sgn * a * h))));
        candidates.insert(candidates.end(), vs.begin(), vs.end());

        bool found = false;
        T best = -std::numeric_limits<T>::max();
        for (auto& c : candidates) {
            T fc_ = ls(c);
            if ((side == 0) != (fc_ < 0))
                continue;
            T margin = std::numeric_limits<T>::max();
            for (auto& f : P.faces) {
                T d = f.normal.dot((f.vertices[0] - c).to_vector());
                if (d < -tol)
                    margin = -1;
                if (f.loc == where || f.loc == location::ON_INTERFACE)
                    margin = std::min(margin, d / h);
            }
            if (margin < -1e-10)
                continue;
            for (auto& x : P.interface) {
                vec3<T> g = ls.gradient(x);
                T d = sgn * (x - c).to_vector().dot(g / g.norm());
                margin = std::min(margin, d / h);
            }
            if (margin >= -1e-10 && margin > best) {
                best = margin;
                P.centers[side] = c;
                found = true;
            }
        }
        if (!found)
            return false;
    }

    return true;
}

/* Parts of a convex polyhedron on both sides of the plane n.x = d (faces
 * only): the first one is on the side n.x <= d */
template<typename T>
std::array<cut_polyhedron<T>, 2>
clip_polyhedron(const cut_polyhedron<T>& P, const vec3<T>& n, T d)
{
    std::array<cut_polyhedron<T>, 2> ret;
    /* point of the edge [a, b] on the plane (from the smaller end) */
    auto plane_point = [&](point<T,3> a, point<T,3> b) {
        if (point_less(b, a))
            std::swap(a, b);
        T fa = n.dot(a.to_vector()) - d, fb = n.dot(b.to_vector()) - d;
        return a + (b - a) * (fa / (fa - fb));
    };
    for (size_t side = 0; side < 2; side++) {
        auto in = [&](const point<T,3>& x) {
            T f = n.dot(x.to_vector()) - d;
            return (side == 0) ? f <= 0 : f >= 0;
        };
        std::vector<point<T,3>> cap;
        auto add_cap = [&](const point<T,3>& x) {
            for (auto& c : cap)
                if (same_point(c, x))
                    return;
            cap.push_back(x);
        };
        for (auto& f : P.faces) {
            std::vector<point<T,3>> v;
            size_t nv = f.vertices.size();
            for (size_t i = 0; i < nv; i++) {
                const auto& a = f.vertices[i];
                const auto& b = f.vertices[(i+1) % nv];
                if (in(a)) {
                    v.push_back(a);
                    if (n.dot(a.to_vector()) - d == 0)
                        add_cap(a);
                }
                if (in(a) != in(b) && (n.dot(a.to_vector()) - d) != 0 && (n.dot(b.to_vector()) - d) != 0) {
                    auto x = plane_point(a, b);
                    v.push_back(x);
                    add_cap(x);
                }
            }
            std::vector<point<T,3>> w;
            for (auto& p : v)
                if (w.empty() || !same_point(w.back(), p))
                    w.push_back(p);
            while (w.size() > 1 && same_point(w.front(), w.back()))
                w.pop_back();
            if (w.size() < 3)
                continue;
            cut_polyhedron_face<T> nf;
            nf.vertices = w;
            nf.normal = f.normal;
            ret[side].faces.push_back(nf);
        }
        if (cap.size() >= 3) {
            point<T,3> c;
            for (auto& p : cap)
                c = c + p;
            c = c / T(cap.size());
            vec3<T> e1 = (cap[0] - c).to_vector();
            e1 -= e1.dot(n) * n;
            e1 /= e1.norm();
            vec3<T> e2 = n.cross(e1);
            std::sort(cap.begin(), cap.end(), [&](const point<T,3>& a, const point<T,3>& b) {
                vec3<T> ra = (a - c).to_vector(), rb = (b - c).to_vector();
                return std::atan2(ra.dot(e2), ra.dot(e1)) < std::atan2(rb.dot(e2), rb.dot(e1));
            });
            cut_polyhedron_face<T> nf;
            nf.vertices = cap;   /* counterclockwise around n */
            nf.normal = n;
            if (side == 1) {
                std::reverse(nf.vertices.begin(), nf.vertices.end());
                nf.normal = -n;
            }
            ret[side].faces.push_back(nf);
        }
    }
    return ret;
}

struct cut_polyhedron_options {
    size_t max_depth = 3;           // splittings of the polyhedra whose parts are not star-shaped
    size_t refine_levels = 8;       // polyline of the faces where the implicit interface fails
    size_t subdivision_levels = 3;  // fallback
};

/* Cut of a convex polyhedron given by its faces (vertices counterclockwise
 * seen from outside, outward normals): location, cut faces, interface patch;
 * the polyhedra where the patch cannot be built are split by the plane
 * through the middle of their interface, orthogonal to its largest extent,
 * up to max_depth times, then integrated by subdivision. */
template<typename T, typename LS>
void
build_cut_polyhedron(cut_polyhedron<T>& P, const LS& ls, const cut_polyhedron_options& opt, size_t depth = 0)
{
    P.children.clear();
    P.loop.clear();
    P.interface.clear();
    P.linear = false;
    P.cut_more_than_once = false;

    std::vector<point<T,3>> all;
    for (auto& f : P.faces)
        all.insert(all.end(), f.vertices.begin(), f.vertices.end());
    P.h = 0;
    for (size_t i = 0; i < all.size(); i++)
        for (size_t j = i+1; j < all.size(); j++)
            P.h = std::max(P.h, (all[i] - all[j]).to_vector().norm());

    std::vector<std::vector<location>> locs(P.faces.size());
    bool any_cut = false;
    for (size_t i = 0; i < P.faces.size(); i++) {
        auto& f = P.faces[i];
        bool has_neg = false, has_pos = false;
        for (auto& v : f.vertices) {
            bool neg = ls(v) < 0;
            locs[i].push_back(neg ? location::IN_NEGATIVE_SIDE : location::IN_POSITIVE_SIDE);
            has_neg = has_neg || neg;
            has_pos = has_pos || !neg;
        }
        f.loc = (has_neg && has_pos) ? location::ON_INTERFACE
              : (has_neg ? location::IN_NEGATIVE_SIDE : location::IN_POSITIVE_SIDE);
        any_cut = any_cut || (f.loc == location::ON_INTERFACE);
    }
    if (!any_cut) {
        P.loc = P.faces.empty() ? location::IN_POSITIVE_SIDE : P.faces[0].loc;
        return;
    }
    P.loc = location::ON_INTERFACE;

    std::map<edge_key<T>, point<T,3>> edge_points;
    bool ok = true;
    for (size_t i = 0; i < P.faces.size() && ok; i++)
        if (P.faces[i].loc == location::ON_INTERFACE)
            ok = make_cut_face(P.faces[i], locs[i], ls, edge_points, opt.refine_levels);
    if (ok && make_interface_patch(P, ls))
        return;

    if (depth >= opt.max_depth || P.faces.size() < 4) {
        P.linear = true;
        return;
    }

    /* splitting plane: through the middle of the interface (samples of the
     * interface curves of the cut faces), orthogonal to its largest extent */
    std::vector<point<T,3>> samples;
    for (auto& f : P.faces) {
        if (f.loc != location::ON_INTERFACE || f.poly.polyline.empty())
            continue;
        planar_level_set<T,3> phi(ls, f.origin, f.U);
        for (auto& p : f.poly.polyline)
            samples.push_back(phi.to_global(p));
    }
    if (samples.size() < 2)
        samples = all;
    point<T,3> c;
    for (auto& p : samples)
        c = c + p;
    c = c / T(samples.size());
    Eigen::Matrix<T,3,3> C = Eigen::Matrix<T,3,3>::Zero();
    for (auto& p : samples) {
        vec3<T> r = (p - c).to_vector();
        C += r * r.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,3,3>> es(C);
    vec3<T> n = es.eigenvectors().col(2);
    T d = n.dot(c.to_vector());

    auto parts = clip_polyhedron(P, n, d);
    if (parts[0].faces.size() < 4 || parts[1].faces.size() < 4) {
        P.linear = true;
        return;
    }
    for (auto& part : parts) {
        build_cut_polyhedron(part, ls, opt, depth + 1);
        P.children.push_back(part);
    }
}

/* Tetrahedra of a convex polyhedron: its centroid and the triangles of its faces */
template<typename T>
std::vector<std::array<point<T,3>, 4>>
polyhedron_tetrahedra(const cut_polyhedron<T>& P)
{
    std::vector<std::array<point<T,3>, 4>> ret;
    point<T,3> c;
    size_t nv = 0;
    for (auto& f : P.faces)
        for (auto& v : f.vertices) {
            c = c + v;
            nv++;
        }
    c = c / T(nv);
    for (auto& f : P.faces) {
        const auto& v = f.vertices;
        if (v.size() == 3) {
            ret.push_back({c, v[0], v[1], v[2]});
            continue;
        }
        point<T,3> of;
        for (auto& p : v)
            of = of + p;
        of = of / T(v.size());
        for (size_t i = 0; i < v.size(); i++)
            ret.push_back({c, of, v[i], v[(i+1) % v.size()]});
    }
    return ret;
}

/* Quadrature of the interface of a cut polyhedron that is not split: the
 * points of the sectors of D lifted onto the interface; weights dA |grad
 * phi| / |grad phi . nu|, normals grad phi / |grad phi|. Returns false if a
 * point cannot be lifted, or the interface is not a graph above D there. */
template<typename T, typename LS>
bool
patch_quadrature(const cut_polyhedron<T>& P, const LS& ls, size_t degree, cut_quadratures<T,3>& out)
{
    using point3 = point<T,3>;
    const auto& o  = P.plane_origin;
    const auto& cD = P.D_center;
    auto lift_plane = [&](const point<T,2>& y) {
        return o + point3(vec3<T>(P.a1 * y.x() + P.a2 * y.y()));
    };

    const auto& rule_xi = gauss_legendre_01<T>((degree + 1) / 2 + 3);
    for (auto& entry : P.loop) {
        face_interface_curve<T> fcrv(ls, P.faces[entry.first]);
        size_t np = fcrv.pieces();
        size_t nt = fcrv.curved() ? (degree + 1) / 2 + 6 : (degree + 1) / 2 + 3;
        for (size_t kp = 0; kp < np; kp++) {
            size_t k = entry.second ? np - 1 - kp : kp;
            for (auto& tw : gauss_legendre_01<T>(nt)) {
                T t = entry.second ? 1.0 - tw.first : tw.first;
                point3 X;
                vec3<T> dX;
                fcrv.eval(k, t, X, dX);
                if (entry.second)
                    dX = -dX;
                vec3<T> r = (X - o).to_vector();
                point<T,2> b(r.dot(P.a1), r.dot(P.a2));
                vec2<T> db(dX.dot(P.a1), dX.dot(P.a2));
                T sb = r.dot(P.nu);
                auto g = b - cD;
                T det = cross(g, db);
                for (auto& xw : rule_xi) {
                    auto y = cD + g * xw.first;
                    T w = tw.second * xw.second * xw.first * det;
                    T s0 = (1.0 - xw.first) * P.center_height + xw.first * sb;
                    auto x = lift_to_interface(ls, lift_plane(y), P.nu, s0, P.h);
                    if (!x)
                        return false;
                    vec3<T> gr = ls.gradient(*x);
                    T gn = gr.dot(P.nu);
                    if (gn <= 0.2 * gr.norm())
                        return false;
                    out.itf.push_back(make_qp(*x, w * gr.norm() / gn));
                    out.normals.push_back(gr / gr.norm());
                }
            }
        }
    }
    return true;
}

/* Quadratures (degree "degree") of the parts of a polyhedron and of its
 * interface by subdivision into tetrahedra */
template<typename T, typename LS>
void
subdivision_quadratures(const cut_polyhedron<T>& P, const LS& ls, size_t degree, size_t levels,
                        cut_quadratures<T,3>& out)
{
    for (auto& tet : polyhedron_tetrahedra(P))
        subdivide_tetrahedron(tet, ls, levels, degree, out);
}

/* Quadrature of the interface of a cut polyhedron (points, weights,
 * normals pointing towards the positive side) */
template<typename T, typename LS>
void
polyhedron_interface_quadrature(const cut_polyhedron<T>& P, const LS& ls, size_t degree,
                                const cut_polyhedron_options& opt, cut_quadratures<T,3>& out)
{
    if (!P.children.empty()) {
        for (auto& c : P.children)
            polyhedron_interface_quadrature(c, ls, degree, opt, out);
        return;
    }
    if (P.loc != location::ON_INTERFACE)
        return;
    cut_quadratures<T,3> q;
    if (P.linear || !patch_quadrature(P, ls, degree, q)) {
        cut_quadratures<T,3> s;
        subdivision_quadratures(P, ls, degree, opt.subdivision_levels, s);
        q.itf = s.itf;
        q.normals = s.normals;
    }
    out.itf.insert(out.itf.end(), q.itf.begin(), q.itf.end());
    out.normals.insert(out.normals.end(), q.normals.begin(), q.normals.end());
}

/* Quadrature of degree "degree" of the part "where" of a polyhedron: cones
 * from its center over its flat faces and the interface (positive weights,
 * points in the part), the uncut parts by cones from their centroid */
template<typename T, typename LS>
void
polyhedron_volume_quadrature(const cut_polyhedron<T>& P, const LS& ls, size_t degree, location where,
                             const cut_polyhedron_options& opt, std::vector<quadrature_point<T,3>>& out)
{
    if (!P.children.empty()) {
        for (auto& c : P.children)
            polyhedron_volume_quadrature(c, ls, degree, where, opt, out);
        return;
    }
    if (P.loc != location::ON_INTERFACE && P.loc != where)
        return;

    bool neg = (where == location::IN_NEGATIVE_SIDE);
    const auto& rule_xi = gauss_legendre_01<T>((degree + 2) / 2 + 1);
    point<T,3> c;
    auto add_cone = [&](const point<T,3>& y, T w, const vec3<T>& n_out) {
        T hgt = (y - c).to_vector().dot(n_out);
        for (auto& xw : rule_xi)
            out.push_back(make_qp(c + (y - c) * xw.first, w * hgt * xw.first * xw.first * xw.second));
    };

    if (P.loc == where) {
        size_t nv = 0;
        for (auto& f : P.faces)
            for (auto& v : f.vertices) {
                c = c + v;
                nv++;
            }
        c = c / T(nv);
        for (auto& f : P.faces) {
            std::vector<quadrature_point<T,3>> fq;
            add_polygon_quadrature(fq, degree, f.vertices);
            for (auto& qp : fq)
                add_cone(qp.point(), qp.weight(), f.normal);
        }
        return;
    }

    cut_quadratures<T,3> itf;
    bool lifted = !P.linear && patch_quadrature(P, ls, std::max<size_t>(degree, 2), itf);
    if (!lifted) {
        cut_quadratures<T,3> s;
        subdivision_quadratures(P, ls, degree, opt.subdivision_levels, s);
        const auto& part = neg ? s.neg : s.pos;
        out.insert(out.end(), part.begin(), part.end());
        return;
    }

    c = P.centers[neg ? 0 : 1];
    for (auto& f : P.faces) {
        if (f.loc == where) {
            std::vector<quadrature_point<T,3>> fq;
            add_polygon_quadrature(fq, degree, f.vertices);
            for (auto& qp : fq)
                add_cone(qp.point(), qp.weight(), f.normal);
        }
        else if (f.loc == location::ON_INTERFACE) {
            planar_level_set<T,3> phi(ls, f.origin, f.U);
            for (auto& qp : side_quadrature(f.poly, phi, degree, where))
                add_cone(phi.to_global(qp.first), qp.second, f.normal);
        }
    }
    for (size_t q = 0; q < itf.itf.size(); q++)
        add_cone(itf.itf[q].point(), itf.itf[q].weight(), neg ? itf.normals[q] : vec3<T>(-itf.normals[q]));
}

/* Sample points of the interface in a polyhedron (bounding boxes, output) */
template<typename T, typename LS>
void
polyhedron_interface_samples(const cut_polyhedron<T>& P, const LS& ls, std::vector<point<T,3>>& out)
{
    if (!P.children.empty()) {
        for (auto& c : P.children)
            polyhedron_interface_samples(c, ls, out);
        return;
    }
    if (P.loc != location::ON_INTERFACE)
        return;
    if (P.linear) {
        cut_quadratures<T,3> s;
        subdivision_quadratures(P, ls, 0, 1, s);
        for (auto& qp : s.itf)
            out.push_back(qp.point());
        return;
    }
    out.insert(out.end(), P.interface.begin(), P.interface.end());
}

/* Statistics of a cut polyhedron: number of leaves split from it, number of
 * leaves integrated by subdivision */
template<typename T>
void
polyhedron_statistics(const cut_polyhedron<T>& P, size_t& split, size_t& linear)
{
    if (!P.children.empty()) {
        split++;
        for (auto& c : P.children)
            polyhedron_statistics(c, split, linear);
        return;
    }
    if (P.loc == location::ON_INTERFACE && P.linear)
        linear++;
}

} // namespace cutgeo
} // namespace disk
