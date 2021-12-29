#![allow(unused)]

use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::consts::FEQ_EPSILON;

/// Info associated with a Triangle shape.
///
/// Fields `p1`, `p2` and `p3` define the points of the `Triangle` in space.
/// These fields are supplied by the entity instantiating the `Triangle`.
///
/// Fields `e1`, `e2` and `normal` are calculated on instantiation. The `e1`
/// and `e2  fields are edges of the triangle, and `normal` is normal to the
/// face of the triangle.
#[derive(Debug, PartialEq)]
pub struct TriangleInfo {
    pub p1: Tuple4D,
    pub p2: Tuple4D,
    pub p3: Tuple4D,

    pub e1: Tuple4D,
    pub e2: Tuple4D,
    pub normal: Tuple4D,
}

impl TriangleInfo {
    pub(crate) fn new(p1: Tuple4D, p2: Tuple4D, p3: Tuple4D) -> TriangleInfo {
        let e1 = p2 - p1;
        let e2 = p3 - p1;
        let normal = e2.cross(&e1).normalize();
        TriangleInfo { p1, p2, p3, e1, e2, normal }
    }
}

#[derive(Debug, PartialEq)]
pub struct SmoothTriangleInfo {
    pub triangle_info: TriangleInfo,

    pub n1: Tuple4D,
    pub n2: Tuple4D,
    pub n3: Tuple4D,
}

impl SmoothTriangleInfo {
    pub(crate) fn new(p1: Tuple4D, p2: Tuple4D, p3: Tuple4D,
        n1: Tuple4D, n2: Tuple4D, n3: Tuple4D) -> SmoothTriangleInfo {

        SmoothTriangleInfo {
            triangle_info: TriangleInfo::new(p1, p2, p3),
            n1, n2, n3
        }
    }
}

/// A bounding box for a shape.
///
/// Helps prevent unnecessary intersection calculations by immediately telling
/// a ray that it lies outside of a valid intersection region.
///
/// For example, imagine a scene with a single sphere (nothing else). Without 
/// bounding boxes, rays would attempt to intersect every pixel in the scene,
/// even though only one sphere is provided. With bounding boxes, rays would
/// terminate intersection calculations immediately if they are unable to find a
/// bounding box.
#[derive(Clone, Debug, PartialEq)]
pub struct Bounds {
    pub minimum: Tuple4D,
    pub maximum: Tuple4D,
}

impl Bounds {
    pub(crate) fn new(min_x: f64, min_y: f64, min_z: f64,
        max_x: f64, max_y: f64, max_z: f64) -> Bounds {
        Bounds {
            minimum: Tuple4D::point(min_x, min_y, min_z),
            maximum: Tuple4D::point(max_x, max_y, max_z)
        }
    }

    /// Creates bounds which encapsulate nothing.
    pub(crate) fn empty() -> Bounds {
        Bounds {
            minimum: Tuple4D::point(0.0, 0.0, 0.0),
            maximum: Tuple4D::point(0.0, 0.0, 0.0)
        }
    }

    /// Creates bounds which encapsulate everything.
    pub(crate) fn infinite() -> Bounds {
        let inf = std::f64::INFINITY;
        Bounds {
            minimum: Tuple4D::point(-inf, -inf, -inf),
            maximum: Tuple4D::point( inf,  inf,  inf)
        }
    }

    /// Transforms a bounding box.
    ///
    /// Certain transformations (e.g. rotations) can cause the bounding box to
    /// expand, as we want the box to be aligned to axis.
    ///
    /// As a result, this function not only transforms the corners of the
    /// existing bounds (represented by `self`), but it also selects the minimum
    /// and maximum components of the transformed corners, producing new bounds.
    pub(crate) fn transform(&self, m: &Matrix4D) -> Bounds {
        let corners = [
            *m * self.minimum,
            *m * Tuple4D::point(self.minimum.x, self.minimum.y, self.maximum.z),
            *m * Tuple4D::point(self.minimum.x, self.maximum.y, self.minimum.z),
            *m * Tuple4D::point(self.minimum.x, self.maximum.y, self.maximum.z),
            *m * Tuple4D::point(self.maximum.x, self.minimum.y, self.minimum.z),
            *m * Tuple4D::point(self.maximum.x, self.minimum.y, self.maximum.z),
            *m * Tuple4D::point(self.maximum.x, self.maximum.y, self.minimum.z),
            *m * self.maximum,
        ];

        let mut min_x = std::f64::INFINITY;
        let mut min_y = std::f64::INFINITY;
        let mut min_z = std::f64::INFINITY;
        let mut max_x = -1.0 * std::f64::INFINITY;
        let mut max_y = -1.0 * std::f64::INFINITY;
        let mut max_z = -1.0 * std::f64::INFINITY;

        for corner in corners.iter() {
            min_x = min_x.min(corner.x);
            min_y = min_y.min(corner.y);
            min_z = min_z.min(corner.z);
            max_x = max_x.max(corner.x);
            max_y = max_y.max(corner.y);
            max_z = max_z.max(corner.z);
        }

        Bounds {
            minimum: Tuple4D::point(min_x, min_y, min_z),
            maximum: Tuple4D::point(max_x, max_y, max_z)
        }
    }

    /// Gets the min and max intersection offsets along an axis of a bound.
    ///
    /// No particular axis is specified; this function takes a component from
    /// the `origin` and `direction` fields of a `Ray4D` (e.g. `origin.x` and
    /// `direction.x`) and returns where the `Ray4D` intersects planes on a cube
    /// for a single axis.
    ///
    /// The smaller `t` is first in the tuple, the larger `t` is second.
    pub(crate) fn check_axis(min: f64, max: f64, origin: f64, direction: f64)
        -> (f64, f64) {
        let tmin_numerator = min - origin;
        let tmax_numerator = max - origin;

        // Make sure that the direction is non-zero. If it is, assign INFINITY.
        let (tmin, tmax) = if direction.abs() >= FEQ_EPSILON {
            (tmin_numerator / direction, tmax_numerator / direction) 
        } else {
            (tmin_numerator * std::f64::INFINITY,
             tmax_numerator * std::f64::INFINITY)
        };

        // If tmin is actually greater than tmax, return tmax first.
        if tmin > tmax {
            (tmax, tmin)
        } else {
            (tmin, tmax)
        }
    }
}

/// The default bounding box should encompass everything.
impl Default for Bounds {
    fn default() -> Bounds {
        let neg_inf = -1.0 * std::f64::INFINITY;
        let pos_inf = std::f64::INFINITY;

        Bounds {
            minimum: Tuple4D::point(neg_inf, neg_inf, neg_inf),
            maximum: Tuple4D::point(pos_inf, pos_inf, pos_inf)
        }
    }
}
