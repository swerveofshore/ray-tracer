use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;

/// A ray in 3D space.
///
/// A ray has an origin (a point with `w == 1.0`) and a direction (a vector with
/// `w == 0.0`). In this crate, rays are primarily cast from the `Camera` onto
/// `Shape`s in the `World` to determine the color of pixels on a `Canvas`.
///
/// See the module documentation for `camera`, `shape`, and `world` for more
/// information.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Ray4D {
    pub origin: Tuple4D,
    pub direction: Tuple4D,
}

impl Ray4D {
    /// Create a new ray with origin and direction.
    pub fn new(mut origin: Tuple4D, mut direction: Tuple4D) -> Ray4D {
        if !origin.is_point() {
            origin.w = 1.0;
        }

        if !direction.is_vector() {
            direction.w = 0.0;
        }

        Ray4D { origin, direction }
    }

    /// Find the position of a ray at time `t`.
    ///
    /// Returns a point using the following formula:
    ///
    /// ```text
    /// point = ray.origin + (t * ray.direction)
    /// ```
    ///
    /// In other words, this function returns a point offset from the ray's
    /// origin by `t` units of the ray's direction vector.
    pub fn position(&self, t: f64) -> Tuple4D {
        self.origin + (t * self.direction)
    }

    /// Transforms a ray, producing a new ray.
    ///
    /// Applies transformation `m` to both the `origin` and `direction`
    /// components of the ray.
    pub fn transform(&self, m: Matrix4D) -> Ray4D {
        Ray4D {
            origin: m * self.origin,
            direction: m * self.direction,
        }
    }
}

#[test]
fn ray_position() {
    let r = Ray4D::new(
                Tuple4D::point(2.0, 3.0, 4.0),
                Tuple4D::vector(1.0, 0.0, 0.0)
            );

    assert_eq!(r.position(0.0), Tuple4D::point(2.0, 3.0, 4.0));
    assert_eq!(r.position(1.0), Tuple4D::point(3.0, 3.0, 4.0));
    assert_eq!(r.position(-1.0), Tuple4D::point(1.0, 3.0, 4.0));
    assert_eq!(r.position(2.5), Tuple4D::point(4.5, 3.0, 4.0));
}

#[test]
fn ray_translation() {
    let r = Ray4D::new(
                Tuple4D::point(1.0, 2.0, 3.0),
                Tuple4D::vector(0.0, 1.0, 0.0)
            );
    let m = Matrix4D::translation(3.0, 4.0, 5.0);
    let t = r.transform(m);

    assert_eq!(t.origin, Tuple4D::point(4.0, 6.0, 8.0));
    assert_eq!(t.direction, Tuple4D::vector(0.0, 1.0, 0.0));
}

#[test]
fn ray_scaling() {
    let r = Ray4D::new(
                Tuple4D::point(1.0, 2.0, 3.0),
                Tuple4D::vector(0.0, 1.0, 0.0)
            );
    let m = Matrix4D::scaling(2.0, 3.0, 4.0);
    let t = r.transform(m);

    assert_eq!(t.origin, Tuple4D::point(2.0, 6.0, 12.0));
    assert_eq!(t.direction, Tuple4D::vector(0.0, 3.0, 0.0));
}
