use std::ops::{ Add, Sub, Neg, Mul };
use crate::feq;

/// A 4D tuple.
///
/// Contains four floating-point values specifying either a point or direction
/// in 3D space.
///
/// If `w` is near `1.0`, a `Tuple4D` is a point. If `w` is near `0.0`, a
/// `Tuple4D` is a vector.
///
/// # Examples
///
/// ```
/// # #![allow(unused)]
/// let point = Tuple4D::point(1.0, 2.0, 3.0);
/// let vector = Tuple4D::vector(4.0, 5.0, 6.0);
/// ```
#[derive(Debug, Default, Copy, Clone, PartialOrd)]
pub struct Tuple4D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64
}

/// Partial equality for a `Tuple4D`.
///
/// Each component of the `Tuple4D` is compared with the function `feq`
/// (floating point equality). This means that equality is calculated
/// approximately; see the `feq` documentation for more info.
impl PartialEq for Tuple4D {
    fn eq(&self, other: &Tuple4D) -> bool {
        feq(self.x, other.x) &&
            feq(self.y, other.y) &&
            feq(self.z, other.z) &&
            feq(self.w, other.w)
    }
}

/// Conversion from a vector to a `Tuple4D`.
///
/// Takes the first `n` elements of a vector `v` and uses them as components of
/// a `Tuple4D`. If `n == 4`, then the `Tuple4D` will be fully populated, with
/// `x == v[0]` through `w == v[3]`.
///
/// If `n > 4`, elements after `v[3]` will be ignored. If `n < 4`, the `Tuple4D`
/// component defaults are used in place.
impl From<&Vec<f64>> for Tuple4D {
    fn from(v: &Vec<f64>) -> Tuple4D {
        match v.len() {
            0 => Default::default(),
            1 => Tuple4D { x: v[0], ..Default::default() },
            2 => Tuple4D { x: v[0], y: v[1], ..Default::default() },
            3 => Tuple4D { x: v[0], y: v[1], z: v[2], ..Default::default() },
            _ => Tuple4D { x: v[0], y: v[1], z: v[2], w: v[3] },
        }
    }
}

impl Tuple4D {
    /// Creates a tuple with supplied parameters.
    pub fn tuple(x: f64, y: f64, z: f64, w: f64) -> Tuple4D {
        Tuple4D { x, y, z, w }
    }

    /// Creates a point with supplied parameters.
    ///
    /// For all points, `w` is near `1.0`.
    pub fn point(x: f64, y: f64, z: f64) -> Tuple4D {
        Tuple4D { x, y, z, w: 1.0 }
    }

    /// Creates a vector with supplied parameters.
    ///
    /// For all vectors, `w` is near `0.0`.
    pub fn vector(x: f64, y: f64, z: f64) -> Tuple4D {
        Tuple4D { x, y, z, w: 0.0 }
    }

    /// Returns whether a `Tuple4D` is a point.
    pub fn is_point(&self) -> bool {
        feq(self.w, 1.0)
    }

    /// Returns whether a `Tuple4D` is a vector.
    pub fn is_vector(&self) -> bool {
        feq(self.w, 0.0)
    }

    /// Returns the magnitude of a `Tuple4D`.
    ///
    /// This includes the `w` component; the following formula is used:
    ///
    /// ```latex
    /// \sqrt{x^2 + y^2 + z^2 + w^2}
    /// ```
    pub fn magnitude(&self) -> f64 {
        f64::sqrt(
            self.x.powi(2) 
            + self.y.powi(2)
            + self.z.powi(2)
            + self.w.powi(2)
        )
    }

    /// Returns a normalized `Tuple4D`.
    ///
    /// Does not modify the underlying `Tuple4D`. All components are modified,
    /// including the `w` component. Effectively, every component is divided by
    /// the magnitude of the existing `Tuple4D`, causing the new magnitude to
    /// be `1.0`.
    pub fn normalize(&self) -> Tuple4D {
        let mag = self.magnitude();

        Tuple4D {
            x: self.x * (1.0 / mag),
            y: self.y * (1.0 / mag),
            z: self.z * (1.0 / mag),
            w: self.w * (1.0 / mag),
        }
    }

    /// Computes the dot product of two `Tuple4D`s.
    ///
    /// The following formula is used (let `self` be `lhs`, and the other
    /// `Tuple4D` be `rhs`):
    ///
    /// ```text
    /// (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w)
    /// ```
    pub fn dot(&self, other: &Tuple4D) -> f64 {
        self.x * other.x
            + self.y * other.y
            + self.z * other.z
            + self.w * other.w
    }

    /// Computes the cross product of two `Tuple4D`s.
    ///
    /// The cross product returns a vector which is orthogonal to both
    /// `Tuple4D`s. "Orthogonal" means that the dot product of the two
    /// `Tuple4D`s  is near `0.0`.
    ///
    /// Note that `w` is not accounted for in this calculation; since the cross
    /// product produces a *vector*, the resultant `w` is `0.0`.
    pub fn cross(&self, other: &Tuple4D) -> Tuple4D {
        Tuple4D {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
            w: 0.0
        }
    }

    /// Reflects a vector across a normal.
    pub fn reflect(&self, normal: &Tuple4D) -> Tuple4D {
        *self - (*normal * 2.0 * self.dot(normal))
    }
}

/// Adds two `Tuple4D`s together.
///
/// Produces a new `Tuple4D` where each component is the sum of the
/// corresponding components of each operand.
impl Add for Tuple4D {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w
        }
    }
}

/// Subtracts one `Tuple4D` from another.
///
/// Produces a new `Tuple4D` where each component is the difference of the
/// corresponding components of each operand.
impl Sub for Tuple4D {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w
        }
    }
}

/// Negates a `Tuple4D`.
///
/// Each component is negated.
impl Neg for Tuple4D {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w
        }
    }
}

/// Implements scalar right-multiplication for a 4D tuple.
///
/// Effectively, this looks like the following:
///
/// ```
/// use ray_tracer_challenge::tuple::Tuple4D;
///
/// let t = Tuple4D::tuple(1.0, 2.0, 3.0, 4.0);
/// let s = 5.0;
///
/// // (notice how the scalar is on the right)
/// assert_eq!(t * s, Tuple4D::tuple(5.0, 10.0, 15.0, 20.0));
/// ```
impl Mul<f64> for Tuple4D {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
            w: self.w * other
        }
    }
}

/// Implements scalar left-multiplication for a 4D tuple.
///
/// Effectively, this looks like the following:
///
/// ```rust
/// use ray_tracer_challenge::tuple::Tuple4D;
///
/// let t = Tuple4D::tuple(1.0, 2.0, 3.0, 4.0);
/// let s = 5.0;
///
/// // (notice how the scalar is on the left)
/// assert_eq!(s * t, Tuple4D::tuple(5.0, 10.0, 15.0, 20.0));
/// ```
impl Mul<Tuple4D> for f64 {
    type Output = Tuple4D;

    fn mul(self, other: Tuple4D) -> Tuple4D {
        Tuple4D {
            x: self * other.x,
            y: self * other.y,
            z: self * other.z,
            w: self * other.w
        }
    }
}

#[test]
fn add_tuples() {
    let a1 = Tuple4D::tuple(3.0, -2.0, 5.0, 1.0);
    let a2 = Tuple4D::tuple(-2.0, 3.0, 1.0, 0.0);
    
    assert_eq!(a1 + a2, Tuple4D::tuple(1.0, 1.0, 6.0, 1.0));
}

#[test]
fn subtract_points() {
    let p1 = Tuple4D::point(3.0, 2.0, 1.0);
    let p2 = Tuple4D::point(5.0, 6.0, 7.0);

    assert_eq!(p1 - p2, Tuple4D::vector(-2.0, -4.0, -6.0));
}

#[test]
fn subtract_vector_from_point() {
    let p = Tuple4D::point(3.0, 2.0, 1.0);
    let v = Tuple4D::vector(5.0, 6.0, 7.0);

    assert_eq!(p - v, Tuple4D::point(-2.0, -4.0, -6.0));
}

#[test]
fn subtract_vectors() {
    let p1 = Tuple4D::vector(3.0, 2.0, 1.0);
    let p2 = Tuple4D::vector(5.0, 6.0, 7.0);

    assert_eq!(p1 - p2, Tuple4D::vector(-2.0, -4.0, -6.0));
}

#[test]
fn negate_tuple() {
    let a = Tuple4D::tuple(1.0, -2.0, 3.0, -4.0);

    assert_eq!(-a, Tuple4D::tuple(-1.0, 2.0, -3.0, 4.0));
}

#[test]
fn multiply_by_scalar() {
    let a = Tuple4D::tuple(1.0, -2.0, 3.0, -4.0);

    assert_eq!(a * 3.5, Tuple4D::tuple(3.5, -7.0, 10.5, -14.0));
}

#[test]
fn multiply_by_fraction() {
    let a = Tuple4D::tuple(1.0, -2.0, 3.0, -4.0);

    assert_eq!(a * 0.5, Tuple4D::tuple(0.5, -1.0, 1.5, -2.0));
}

#[test]
fn magnitude_of_positive_vector() {
    let v = Tuple4D::vector(1.0, 2.0, 3.0);

    assert_eq!(v.magnitude(), f64::sqrt(14.0));
}

#[test]
fn magnitude_of_negative_vector() {
    let v = Tuple4D::vector(-1.0, -2.0, -3.0);

    assert_eq!(v.magnitude(), f64::sqrt(14.0));
}

#[test]
fn normalize_whole_sum() {
    let v = Tuple4D::vector(4.0, 0.0, 0.0);

    assert_eq!(v.normalize(), Tuple4D::vector(1.0, 0.0, 0.0));
}

#[test]
fn normalize_fractional_sum() {
    let v = Tuple4D::vector(1.0, 2.0, 3.0);
    let e = Tuple4D::vector(
        1.0 / f64::sqrt(14.0),
        2.0 / f64::sqrt(14.0),
        3.0 / f64::sqrt(14.0)
    );

    assert_eq!(v.normalize(), e);
}

#[test]
fn dot_product_of_vectors() {
    let a = Tuple4D::vector(1.0, 2.0, 3.0);
    let b = Tuple4D::vector(2.0, 3.0, 4.0);

    assert_eq!(a.dot(&b), 20.0);
}

#[test]
fn cross_product_of_vectors() {
    let a = Tuple4D::vector(1.0, 2.0, 3.0);
    let b = Tuple4D::vector(2.0, 3.0, 4.0);

    let c = Tuple4D::vector(-1.0, 2.0, -1.0);
    let d = Tuple4D::vector(1.0, -2.0, 1.0);

    assert_eq!(a.cross(&b), c);
    assert_eq!(b.cross(&a), d);
}

#[test]
fn reflect_vector_at_45_degrees() {
    let v = Tuple4D::vector(1.0, -1.0, 0.0);
    let n = Tuple4D::vector(0.0, 1.0, 0.0);
    let r = v.reflect(&n);

    assert_eq!(r, Tuple4D::vector(1.0, 1.0, 0.0));
}
