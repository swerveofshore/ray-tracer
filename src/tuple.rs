use std::ops::{ Add, Sub, Neg, Mul };

use crate::feq;

#[allow(unused)]
#[derive(Debug, Default, Copy, Clone, PartialOrd)]
pub struct Tuple3D {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl PartialEq for Tuple3D {
    fn eq(&self, other: &Tuple3D) -> bool {
        feq(self.x, other.x) &&
            feq(self.y, other.y) &&
            feq(self.z, other.z)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialOrd)]
pub struct Tuple4D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64
}

impl PartialEq for Tuple4D {
    fn eq(&self, other: &Tuple4D) -> bool {
        feq(self.x, other.x) &&
            feq(self.y, other.y) &&
            feq(self.z, other.z) &&
            feq(self.w, other.w)
    }
}

impl Tuple4D {
    pub fn tuple(x: f64, y: f64, z: f64, w: f64) -> Tuple4D {
        Tuple4D { x, y, z, w }
    }

    pub fn point(x: f64, y: f64, z: f64) -> Tuple4D {
        Tuple4D { x, y, z, w: 1.0 }
    }

    pub fn vector(x: f64, y: f64, z: f64) -> Tuple4D {
        Tuple4D { x, y, z, w: 0.0 }
    }

    pub fn is_point(&self) -> bool {
        self.w == 1.0
    }

    pub fn is_vector(&self) -> bool {
        self.w == 0.0
    }

    pub fn magnitude(&self) -> f64 {
        f64::sqrt(
            self.x.powi(2) 
            + self.y.powi(2)
            + self.z.powi(2)
            + self.w.powi(2)
        )
    }

    pub fn normalize(&self) -> Tuple4D {
        let mag = self.magnitude();

        Tuple4D {
            x: self.x * (1.0 / mag),
            y: self.y * (1.0 / mag),
            z: self.z * (1.0 / mag),
            w: self.w * (1.0 / mag),
        }
    }

    pub fn dot(&self, other: &Tuple4D) -> f64 {
        self.x * other.x
            + self.y * other.y
            + self.z * other.z
            + self.w * other.w
    }

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

/* Tests */

#[test]
fn add_tuples() {
    let a1 = Tuple4D::tuple(3.0, -2.0, 5.0, 1.0);
    let a2 = Tuple4D::tuple(-2.0, 3.0, 1.0, 0.0);
    
    assert_eq!(a1 + a2, Tuple4D::tuple(1.0, 1.0, 6.0, 1.0));
}

#[test]
fn sub_points() {
    let p1 = Tuple4D::point(3.0, 2.0, 1.0);
    let p2 = Tuple4D::point(5.0, 6.0, 7.0);

    assert_eq!(p1 - p2, Tuple4D::vector(-2.0, -4.0, -6.0));
}

#[test]
fn sub_vector_from_point() {
    let p = Tuple4D::point(3.0, 2.0, 1.0);
    let v = Tuple4D::vector(5.0, 6.0, 7.0);

    assert_eq!(p - v, Tuple4D::point(-2.0, -4.0, -6.0));
}

#[test]
fn sub_vectors() {
    let p1 = Tuple4D::vector(3.0, 2.0, 1.0);
    let p2 = Tuple4D::vector(5.0, 6.0, 7.0);

    assert_eq!(p1 - p2, Tuple4D::vector(-2.0, -4.0, -6.0));
}

#[test]
fn neg_tuple() {
    let a = Tuple4D::tuple(1.0, -2.0, 3.0, -4.0);

    assert_eq!(-a, Tuple4D::tuple(-1.0, 2.0, -3.0, 4.0));
}

#[test]
fn mul_scalar() {
    let a = Tuple4D::tuple(1.0, -2.0, 3.0, -4.0);

    assert_eq!(a * 3.5, Tuple4D::tuple(3.5, -7.0, 10.5, -14.0));
}

#[test]
fn mul_fraction() {
    let a = Tuple4D::tuple(1.0, -2.0, 3.0, -4.0);

    assert_eq!(a * 0.5, Tuple4D::tuple(0.5, -1.0, 1.5, -2.0));
}

#[test]
fn magnitude_pos() {
    let v = Tuple4D::vector(1.0, 2.0, 3.0);

    assert_eq!(v.magnitude(), f64::sqrt(14.0));
}

#[test]
fn magnitude_neg() {
    let v = Tuple4D::vector(-1.0, -2.0, -3.0);

    assert_eq!(v.magnitude(), f64::sqrt(14.0));
}

#[test]
fn normalize_clean() {
    let v = Tuple4D::vector(4.0, 0.0, 0.0);

    assert_eq!(v.normalize(), Tuple4D::vector(1.0, 0.0, 0.0));
}

#[test]
fn normalize_dirty() {
    let v = Tuple4D::vector(1.0, 2.0, 3.0);
    let e = Tuple4D::vector(
        1.0 / f64::sqrt(14.0),
        2.0 / f64::sqrt(14.0),
        3.0 / f64::sqrt(14.0)
    );

    assert_eq!(v.normalize(), e);
}

#[test]
fn dot_vectors() {
    let a = Tuple4D::vector(1.0, 2.0, 3.0);
    let b = Tuple4D::vector(2.0, 3.0, 4.0);

    assert_eq!(a.dot(&b), 20.0);
}

#[test]
fn cross_vectors() {
    let a = Tuple4D::vector(1.0, 2.0, 3.0);
    let b = Tuple4D::vector(2.0, 3.0, 4.0);

    let c = Tuple4D::vector(-1.0, 2.0, -1.0);
    let d = Tuple4D::vector(1.0, -2.0, 1.0);

    assert_eq!(a.cross(&b), c);
    assert_eq!(b.cross(&a), d);
}

#[test]
fn reflect_45() {
    let v = Tuple4D::vector(1.0, -1.0, 0.0);
    let n = Tuple4D::vector(0.0, 1.0, 0.0);
    let r = v.reflect(&n);

    assert_eq!(r, Tuple4D::vector(1.0, 1.0, 0.0));
}
