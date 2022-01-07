use std::fmt;
use std::ops::{ Index, IndexMut, Mul };
use std::convert::From;

use crate::feq;
use crate::tuple::*;

/// A 2x2 matrix.
#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd)]
struct Matrix2D {
    data: [f64; 4],
}

impl From<[f64; 4]> for Matrix2D {
    fn from(data: [f64; 4]) -> Matrix2D {
        Matrix2D { data } 
    }
}

impl Index<(usize, usize)> for Matrix2D {
    type Output = f64;

    fn index<'a>(&'a self, index: (usize, usize)) -> &'a f64 {
        &self.data[(index.0 * 2) + index.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix2D {
    fn index_mut<'a>(&'a mut self, index: (usize, usize)) -> &'a mut f64 {
        &mut self.data[(index.0 * 2) + index.1]
    }
}

/// A 3x3 matrix.
#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd)]
struct Matrix3D {
    data: [f64; 9],
}

impl From<[f64; 9]> for Matrix3D {
    fn from(data: [f64; 9]) -> Matrix3D {
        Matrix3D { data } 
    }
}

impl Index<(usize, usize)> for Matrix3D {
    type Output = f64;

    fn index<'a>(&'a self, index: (usize, usize)) -> &'a f64 {
        &self.data[(index.0 * 3) + index.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix3D {
    fn index_mut<'a>(&'a mut self, index: (usize, usize)) -> &'a mut f64 {
        &mut self.data[(index.0 * 3) + index.1]
    }
}

/// A 4x4 matrix.
///
/// These matrices are used almost universally in the ray tracer logic.
/// Basically, these matrices encode transformations in 3D space, transforming 
/// both vectors and points (`w` components of `0.0` and `1.0`, respectively).
///
/// Note that smaller matrix types like `Matrix2D` and `Matrix3D` are defined
/// in the same module as `Matrix4D`. However, types `Matrix2D` and `Matrix3D`
/// are seldom used, so they are private. This may change in later versions.
///
/// For methods which modify matrices, they are typically provided in pairs;
/// one which modifies the matrix in-place, and one which returns a new matrix.
/// For example, `transpose` and `transposition`. Method `transpose` takes a
/// `Matrix4D`, and turns said matrix into its own transpose. Method
/// `transposition` produces a new matrix, the transpose of the original matrix.
/// 
/// See the individual method documentation comments for more information.
///
/// # Examples
///
/// Creating an identity matrix:
///
/// ```
/// # #![allow(unused)]
/// # use ray_tracer_challenge::matrix::Matrix4D;
/// let mat = Matrix4D::identity();
/// assert_eq!(mat.determinant(), 1.0);
/// ```
///
/// Calculating a view transformation (for cameras, etc.):
///
/// ```
/// # #![allow(unused)]
/// # use ray_tracer_challenge::tuple::Tuple4D;
/// # use ray_tracer_challenge::matrix::Matrix4D;
/// let from = Tuple4D::point(0.0, 0.0, 0.0);
/// let to = Tuple4D::point(0.0, 0.0, 5.0);
/// let up = Tuple4D::vector(0.0, 1.0, 0.0);
/// let view = Matrix4D::view_transform(from, to, up);
/// ```
#[derive(Copy, Clone, Debug, Default, PartialOrd)]
pub struct Matrix4D {
    data: [f64; 16],
}

impl Matrix2D {
    /// Creates a new `Matrix2D`. All elements are initialized to `0.0`.
    #[allow(unused)]
    fn new() -> Matrix2D {
        Matrix2D { data: [0.0; 4] }
    }

    /// Instantiates a 2x2 identity matrix.
    #[allow(unused)]
    fn identity() -> Matrix2D {
        let mut buf = [0.0; 4];
        buf[0] = 1.0; buf[3] = 1.0;

        Matrix2D { data: buf }
    }

    /// Calculates the determinant of a `Matrix2D`.
    fn determinant(&self) -> f64 {
        self[(0, 0)] * self[(1, 1)] - self[(0, 1)] * self[(1, 0)]
    }
}

impl Matrix3D {
    /// Creates a new `Matrix3D`. All elements are initialized to `0.0`.
    #[allow(unused)]
    fn new() -> Matrix3D {
        Matrix3D { data: [0.0; 9] }
    }

    /// Instantiates a 3x3 identity matrix.
    #[allow(unused)]
    fn identity() -> Matrix3D {
        let mut buf = [0.0; 9];
        buf[0] = 1.0; buf[4] = 1.0; buf[8] = 1.0;

        Matrix3D { data: buf }
    }

    /// Returns the submatrix of a `Matrix3D`.
    ///
    /// A submatrix can be thought of as a matrix which "eliminates" a row and
    /// column of a larger matrix. For example, given the following 3x3 matrix:
    ///
    /// ```text
    /// [
    ///     1.0, 0.0, 2.0,
    ///     3.0, 1.0, 0.0,
    ///     1.0, 1.0, 1.0
    /// ]
    /// ```
    ///
    /// The corresponding submatrix for `row == 1`, `col == 2` (assuming zero
    /// index), would be a 2x2 matrix:
    ///
    /// ```text
    /// [
    ///     1.0, 0.0,
    ///     1.0, 1.0
    /// ]
    /// ```
    ///
    /// (Notice how row index `1` and column index `2` of the original 3x3
    /// matrix are removed, yielding the above 2x2 matrix).
    fn submatrix(&self, row: usize, col: usize) -> Matrix2D {
        let mut buf: [f64; 4] = [0.0; 4];
        let mut count = 0;

        for r in 0..3 {
            for c in 0..3 {
                if !(r == row || c == col) {
                    buf[count] = self[(r, c)];
                    count += 1;
                }
            }
        }

        Matrix2D { data: buf }
    }

    /// Returns the minor of a `Matrix3D` at row and column.
    ///
    /// The "minor" is the determinant of the submatrix at `row` and `col`. See
    /// the documentation for `submatrix` for what this means.
    pub fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    /// Returns the cofactor of a `Matrix3D` at row and column.
    ///
    /// The "cofactor" is the minor of a matrix, negated according to the
    /// "cofactor matrix." Basically, if the sum of row and column is even,
    /// the minor remains positive; if the sum is odd, the minor is negated.
    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let m = self.minor(row, col);
        m * if (row + col) % 2 == 0 { 1.0 } else { -1.0 }
    }

    /// Calculates the determinant of a Matrix3D.
    pub fn determinant(&self) -> f64 {
        let mut sum = 0.0;
        for c in 0..3 {
            sum += self[(0, c)] * self.cofactor(0, c);
        }

        sum
    }
}

/// Determines whether two `Matrix4D`s are equal.
///
/// Matrices are compared element-wise. Note that equality is approximate, as
/// `Matrix4D` elements are floating point numbers.
impl PartialEq for Matrix4D {
    fn eq(&self, other: &Matrix4D) -> bool {
        self.data.iter().zip(other.data.iter()).all(|(x, y)| feq(*x, *y))
    }
}

impl Matrix4D {
    /// Creates a new `Matrix4D`. All elements are initialized to `0.0`.
    pub fn new() -> Matrix4D {
        Matrix4D { data: [0.0; 16] }
    }

    /// Instantiates a 4x4 identity matrix.
    pub fn identity() -> Matrix4D {
        let mut buf = [0.0; 16];
        buf[0] = 1.0; buf[5] = 1.0; buf[10] = 1.0; buf[15] = 1.0;

        Matrix4D { data: buf }
    }

    /// Instantiates a 4x4 translation matrix.
    ///
    /// This matrix offsets a point by `x`, `y` and `z`.
    pub fn translation(x: f64, y: f64, z: f64) -> Matrix4D {
        let mut trans = Self::identity();
        trans[(0, 3)] = x;
        trans[(1, 3)] = y;
        trans[(2, 3)] = z;

        trans
    }

    /// Instantiates a 4x4 scaling matrix.
    ///
    /// This matrix scales vectors or points by `x`, `y` and `z` along the X, Y
    /// and Z axes, respectively.
    pub fn scaling(x: f64, y: f64, z: f64) -> Matrix4D {
        let mut scale = Self::identity();
        scale[(0, 0)] = x;
        scale[(1, 1)] = y;
        scale[(2, 2)] = z;

        scale
    }

    /// Instantiates a 4x4 rotation matrix, rotating about the X axis.
    ///
    /// Rotations occur clockwise. Assumes that parameter `r` is in radians.
    ///
    /// # Examples
    ///
    /// Create a matrix to rotate a point 90 degrees about the X axis:
    ///
    /// ```
    /// # #![allow(unused)]
    /// # use ray_tracer_challenge::tuple::Tuple4D;
    /// # use ray_tracer_challenge::matrix::Matrix4D;
    /// let point = Tuple4D::point(0.0, 1.0, 0.0);
    /// let m = Matrix4D::rotation_x(std::f64::consts::PI / 2.0);
    /// assert_eq!(m * point, Tuple4D::point(0.0, 0.0, 1.0));
    /// ```
    pub fn rotation_x(r: f64) -> Matrix4D {
        let mut rotate = Self::identity();
        rotate[(1, 1)] =  r.cos();
        rotate[(1, 2)] = -r.sin();
        rotate[(2, 1)] =  r.sin();
        rotate[(2, 2)] =  r.cos();

        rotate
    }

    /// Instantiates a 4x4 rotation matrix, rotating about the Y axis.
    ///
    /// Rotations occur clockwise. Assumes that parameter `r` is in radians.
    ///
    /// # Examples
    ///
    /// Create a matrix to rotate a point 90 degrees about the Y axis:
    ///
    /// ```
    /// # #![allow(unused)]
    /// # use ray_tracer_challenge::tuple::Tuple4D;
    /// # use ray_tracer_challenge::matrix::Matrix4D;
    /// let point = Tuple4D::point(1.0, 0.0, 0.0);
    /// let m = Matrix4D::rotation_y(std::f64::consts::PI / 2.0);
    /// assert_eq!(m * point, Tuple4D::point(0.0, 0.0, -1.0));
    /// ```
    pub fn rotation_y(r: f64) -> Matrix4D {
        let mut rotate = Self::identity();
        rotate[(0, 0)] =  r.cos();
        rotate[(0, 2)] =  r.sin();
        rotate[(2, 0)] = -r.sin();
        rotate[(2, 2)] =  r.cos();

        rotate
    }

    /// Instantiates a 4x4 rotation matrix, rotating about the Z axis.
    ///
    /// Rotations occur clockwise. Assumes that parameter `r` is in radians.
    ///
    /// # Examples
    ///
    /// Create a matrix to rotate a point 90 degrees about the Z axis:
    ///
    /// ```
    /// # #![allow(unused)]
    /// # use ray_tracer_challenge::tuple::Tuple4D;
    /// # use ray_tracer_challenge::matrix::Matrix4D;
    /// let point = Tuple4D::point(0.0, 1.0, 0.0);
    /// let m = Matrix4D::rotation_z(std::f64::consts::PI / 2.0);
    /// assert_eq!(m * point, Tuple4D::point(-1.0, 0.0, 0.0));
    /// ```
    pub fn rotation_z(r: f64) -> Matrix4D {
        let mut rotate = Self::identity();
        rotate[(0, 0)] =  r.cos();
        rotate[(0, 1)] = -r.sin();
        rotate[(1, 0)] =  r.sin();
        rotate[(1, 1)] =  r.cos();

        rotate
    }

    /// Instantiates a 4x4 shearing matrix.
    ///
    /// TODO: Probably ignore this documentation; it needs work. It may be
    /// helpful for understanding the high-level of what's going on, but it
    /// doesn't explain what's happening in a shearing transformation very well.
    ///
    /// ---
    ///
    /// A shearing matrix creates a kind of "slope" along all specified axis.
    /// To elaborate, a shear takes something like this (imagine a square):
    ///
    /// ```text
    /// ------
    /// |    |
    /// |    |
    /// ------
    /// ```
    ///
    /// And turns it into this:
    ///
    /// ```text
    /// ------
    /// \     \
    ///  \     \
    ///   ------
    /// ```
    ///
    /// Basically, some axis is manipulated to make the shape look literally
    /// "sheared," or distorted in some direction.
    ///
    /// The above example alters the change in `x` with respect to `y` for all
    /// points on the square. This is the only possible transformation in 2D,
    /// really.
    ///
    /// In 3D, more axis can change with respect to one another; these changes
    /// are provided as arguments to this function. For example, if we wanted
    /// to alter the change in `x` with respect to `y` in 3D, we would specify
    /// the `xy` parameter as some number above or below `1.0`.
    ///
    /// # Examples
    ///
    /// Shearing points along the `xy` slope (similar to the above example):
    ///
    /// ```
    /// # use ray_tracer_challenge::tuple::Tuple4D;
    /// # use ray_tracer_challenge::matrix::Matrix4D;
    /// let point = Tuple4D::point(1.0, 0.0, 0.0);
    /// let m = Matrix4D::shearing(2.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    /// assert_eq!(m * point, Tuple4D::point(1.0, 1.0, 1.0));
    /// ```
    pub fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64)
        -> Matrix4D {
        let mut shear = Self::identity();
        shear[(0, 1)] = xy;
        shear[(0, 2)] = xz;
        shear[(1, 0)] = yx;
        shear[(1, 2)] = yz;
        shear[(2, 0)] = zx;
        shear[(2, 1)] = zy;

        shear
    }

    /// Generates a view transformation.
    ///
    /// The view transform manipulates the world from the perspective of an eye,
    /// The `from` parameter is where the eye is, the `to` parameter is where
    /// the eye is looking, and the `up` parameter indicates where "up" is in
    /// the world.
    ///
    /// A "default" orientation fixes the eye at the origin, looking at a screen
    /// one unit "deep." The `up` vector points conventionally up, with `y=1`.
    ///
    /// Note that the view transformation moves the *world* with respect to the
    /// eye, not the other way around. 
    pub fn view_transform(from: Tuple4D, to: Tuple4D, up: Tuple4D) -> Matrix4D {
        let forward = (to - from).normalize(); 
        let left = forward.cross(&up.normalize());
        let true_up = left.cross(&forward);

        let mut orientation = Matrix4D::identity();
        orientation[(0, 0)] = left.x;
        orientation[(0, 1)] = left.y;
        orientation[(0, 2)] = left.z;

        orientation[(1, 0)] = true_up.x;
        orientation[(1, 1)] = true_up.y;
        orientation[(1, 2)] = true_up.z;

        orientation[(2, 0)] = -forward.x;
        orientation[(2, 1)] = -forward.y;
        orientation[(2, 2)] = -forward.z;

        orientation * Matrix4D::translation(-from.x, -from.y, -from.z)
    }

    /// Produces the transpose of a matrix in-place.
    ///
    /// The transpose of a matrix is roughly defined by the following formula
    /// (given matrix `A`, create transpose matrix `A^T`):
    ///
    /// ```latex
    /// A^T_{ij} = A_{ji}
    /// ```
    ///
    /// Where subscripts `ij` represent the element of `A` at row `i`, column
    /// `j`.
    pub fn transpose(&mut self) {
        for r in 0..4 {
            for c in 1..(4-r) {
                let tmp = self[(r, c)];
                self[(r, c)] = self[(c, r)];
                self[(c, r)] = tmp;
            }
        }
    }

    /// Produces the transpose of a matrix, returning a new matrix as a result.
    ///
    /// See the documentation on method `transpose` for more information on what
    /// a transpose is.
    pub fn transposition(&self) -> Matrix4D {
        let mut buf = self.clone();

        for r in 0..4 {
            for c in (r+1)..4 {
                let tmp = buf[(r, c)];
                buf[(r, c)] = buf[(c, r)];
                buf[(c, r)] = tmp;
            }
        }

        buf
    }

    /// Returns the submatrix of a `Matrix4D`.
    ///
    /// A submatrix can be thought of as a matrix which "eliminates" a row and
    /// column of a larger matrix. For example, given the following 4x4 matrix:
    ///
    /// ```text
    /// [
    ///     1.0, 0.0, 2.0, 0.0,
    ///     3.0, 1.0, 0.0, 0.0,
    ///     1.0, 1.0, 1.0, 0.0,
    ///     4.0, 0.0, 0.0, 1.0
    /// ]
    /// ```
    ///
    /// The corresponding submatrix for `row == 1`, `col == 2` (assuming zero
    /// index), would be a 3x3 matrix:
    ///
    /// ```text
    /// [
    ///     1.0, 0.0, 0.0,
    ///     1.0, 1.0, 0.0,
    ///     4.0, 0.0, 1.0
    /// ]
    /// ```
    ///
    /// (Notice how row index `1` and column index `2` of the original 4x4
    /// matrix are removed, yielding the above 2x2 matrix).
    fn submatrix(&self, row: usize, col: usize) -> Matrix3D {
        let mut buf: [f64; 9] = [0.0; 9];
        let mut count = 0;

        for r in 0..4 {
            for c in 0..4 {
                if !(r == row || c == col) {
                    buf[count] = self[(r, c)];
                    count += 1;
                }
            }
        }

        Matrix3D { data: buf }
    }

    /// Returns the minor of a `Matrix4D` at row and column.
    ///
    /// The "minor" is the determinant of the submatrix at `row` and `col`. See
    /// the documentation for `submatrix` for what this means.
    pub fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    /// Returns the cofactor of a `Matrix4D` at row and column.
    ///
    /// The "cofactor" is the minor of a matrix, negated according to the
    /// "cofactor matrix." Basically, if the sum of row and column is even,
    /// the minor remains positive; if the sum is odd, the minor is negated.
    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let m = self.minor(row, col);
        m * if (row + col) % 2 == 0 { 1.0 } else { -1.0 }
    }

    /// Calculates the determinant of a `Matrix4D`.
    pub fn determinant(&self) -> f64 {
        let mut sum = 0.0;
        for c in 0..4 {
            sum += self[(0, c)] * self.cofactor(0, c);
        }

        sum
    }

    /// Calculates the inverse of a `Matrix4D`, if it exists.
    ///
    /// If a `Matrix4D` is non-invertible, this function returns `None`. A
    /// `Matrix4D` is determined invertible if its determinant is nonzero.
    pub fn inverse(&self) -> Option<Matrix4D> {
        let det = self.determinant();
        if det == 0.0 {
            return None;
        }

        let mut inv = Matrix4D::new();
        for r in 0..4 {
            for c in 0..4 {
                inv[(c, r)] = self.cofactor(r, c) / det;
            }
        }

        Some(inv)
    }
}

impl From<[f64; 16]> for Matrix4D {
    fn from(data: [f64; 16]) -> Matrix4D {
        Matrix4D { data } 
    }
}

impl Index<(usize, usize)> for Matrix4D {
    type Output = f64;

    fn index<'a>(&'a self, index: (usize, usize)) -> &'a f64 {
        &self.data[(index.0 * 4) + index.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix4D {
    fn index_mut<'a>(&'a mut self, index: (usize, usize)) -> &'a mut f64 {
        &mut self.data[(index.0 * 4) + index.1]
    }
}

/// Multiplication between two matrices.
///
/// Note that matrix multiplication is not commutative; in other words, for
/// matrix `A` and matrix `B`, `A * B` is not necessarily equal to `B * A`.
///
/// Certain transformations like scaling are commutative about themselves
/// (e.g. a scaling matrix multiplied by a scaling matrix), but this is not
/// guaranteed for all matrix products.
///
/// # Examples
///
/// ```
/// # use ray_tracer_challenge::matrix::Matrix4D;
/// let m1 = Matrix4D::scaling(2.0, 3.0, 4.0);
/// let m2 = Matrix4D::scaling(4.0, 3.0, 2.0);
/// assert_eq!(m1 * m2, Matrix4D::scaling(8.0, 9.0, 8.0));
/// ```
impl Mul<Matrix4D> for Matrix4D {
    type Output = Matrix4D;

    fn mul(self, other: Matrix4D) -> Matrix4D {
        let mut res = Matrix4D::new();

        for r in 0..4 {
            for c in 0..4 {
                res[(r, c)] = self[(r, 0)] * other[(0, c)]
                    + self[(r, 1)] * other[(1, c)]
                    + self[(r, 2)] * other[(2, c)]
                    + self[(r, 3)] * other[(3, c)]
            }
        }

        res
    }
}

/// Multiplication between a matrix and a `Tuple4D`.
///
/// Note that `Tuple4D`s are multiplied on the right; this requirement is
/// somewhat arbitrary, but it matches the convention of a 4D vector having
/// 4 rows, 1 column.
///
/// # Examples
///
/// ```
/// # use ray_tracer_challenge::tuple::Tuple4D;
/// # use ray_tracer_challenge::matrix::Matrix4D;
/// let v = Tuple4D::vector(1.0, 4.0, 5.0);
/// let m = Matrix4D::scaling(2.0, 2.0, 2.0);
/// assert_eq!(m * v, Tuple4D::vector(2.0, 8.0, 10.0));
/// ```
impl Mul<Tuple4D> for Matrix4D {
    type Output = Tuple4D;

    fn mul(self, other: Tuple4D) -> Tuple4D {
        let mut buf: [f64; 4] = Default::default();

        for r in 0..4 {
            buf[r] = self[(r, 0)] * other.x
                + self[(r, 1)] * other.y
                + self[(r, 2)] * other.z
                + self[(r, 3)] * other.w;
        }

        Tuple4D { x: buf[0], y: buf[1], z: buf[2], w: buf[3] }
    }
}

impl fmt::Display for Matrix4D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for r in 0..4 {
            write!(f, "|")?;
            for c in 0..4 {
                write!(f, " {} |", self[(r, c)])?;
            }

            // Don't put a newline on the final row (allow the user to do that)
            if r != 3 {
                write!(f, "\n")?;
            }
        }

        Ok(())
    }
}

#[test]
fn identity() {
    let i = Matrix4D::identity();
    let a: Matrix4D = [ 0.0, 1.0,  2.0,  4.0, 
                        1.0, 2.0,  4.0,  8.0,
                        2.0, 4.0,  8.0, 16.0,
                        4.0, 8.0, 16.0, 32.0, ].into();

    assert_eq!(i * a, a);
    assert_eq!(a * i, a);
}

#[test]
fn transpose() {
     let a: Matrix4D = [ 0.0, 9.0, 3.0, 0.0, 
                         9.0, 8.0, 0.0, 8.0,
                         1.0, 8.0, 5.0, 3.0,
                         0.0, 0.0, 5.0, 8.0, ].into();

     let t: Matrix4D = [ 0.0, 9.0, 1.0, 0.0, 
                         9.0, 8.0, 8.0, 0.0,
                         3.0, 0.0, 5.0, 5.0,
                         0.0, 8.0, 3.0, 8.0, ].into();

     assert_eq!(t, a.transposition());
     assert_eq!(t.transposition(), a);
}

#[test]
fn transpose_identity() {
    let i = Matrix4D::identity();
    assert_eq!(i, i.transposition());
}

#[test]
fn mat3_submatrix() {
    let a: Matrix3D = [  1.0, 5.0,  0.0,
                        -3.0, 2.0,  7.0,
                         0.0, 6.0, -3.0, ].into();

    let s: Matrix2D = [ -3.0, 2.0,
                         0.0, 6.0  ].into();

    assert_eq!(a.submatrix(0, 2), s);
}

#[test]
fn mat4_submatrix() {
     let a: Matrix4D = [ -6.0, 1.0,  1.0, 6.0, 
                         -8.0, 5.0,  8.0, 6.0,
                         -1.0, 0.0,  8.0, 2.0,
                         -7.0, 1.0, -1.0, 1.0, ].into();

     let s: Matrix3D = [ -6.0,  1.0, 6.0,
                         -8.0,  8.0, 6.0,
                         -7.0, -1.0, 1.0, ].into();

     assert_eq!(a.submatrix(2, 1), s);
}

#[test]
fn mat3_minor() {
    let a: Matrix3D = [ 3.0,  5.0,  0.0,
                        2.0, -1.0, -7.0,
                        6.0, -1.0,  5.0, ].into();

    assert_eq!(a.minor(1, 0), 25.0);
}

#[test]
fn mat3_cofactor() {
    let a: Matrix3D = [ 3.0,  5.0,  0.0,
                        2.0, -1.0, -7.0,
                        6.0, -1.0,  5.0, ].into();

    assert_eq!(a.minor(0, 0), -12.0);
    assert_eq!(a.cofactor(0, 0), -12.0);
    assert_eq!(a.minor(1, 0), 25.0);
    assert_eq!(a.cofactor(1, 0), -25.0);
}

#[test]
fn mat3_determinant() {
     let a: Matrix3D = [  1.0, 2.0,  6.0,
                         -5.0, 8.0, -4.0,
                          2.0, 6.0,  4.0, ].into();

     assert_eq!(a.cofactor(0, 0), 56.0);
     assert_eq!(a.cofactor(0, 1), 12.0);
     assert_eq!(a.cofactor(0, 2), -46.0);
     assert_eq!(a.determinant(), -196.0);
}

#[test]
fn mat4_determinant() {
     let a: Matrix4D = [ -2.0, -8.0,  3.0,  5.0, 
                         -3.0,  1.0,  7.0,  3.0,
                          1.0,  2.0, -9.0,  6.0,
                         -6.0,  7.0,  7.0, -9.0, ].into();

     assert_eq!(a.cofactor(0, 0), 690.0);
     assert_eq!(a.cofactor(0, 1), 447.0);
     assert_eq!(a.cofactor(0, 2), 210.0);
     assert_eq!(a.cofactor(0, 3), 51.0);
     assert_eq!(a.determinant(), -4071.0);
}

#[test]
fn mat4_inverse() {
     let a: Matrix4D = [  8.0, -5.0,  9.0,  2.0, 
                          7.0,  5.0,  6.0,  1.0,
                         -6.0,  0.0,  9.0,  6.0,
                         -3.0,  0.0, -9.0, -4.0, ].into();

     let i: Matrix4D = [ -0.15385, -0.15385, -0.28205, -0.53846,
                         -0.07692,  0.12308,  0.02564,  0.03077,
                          0.35897,  0.35897,  0.43590,  0.92308,
                         -0.69231, -0.69231, -0.76923, -1.92308, ].into();

     // Fails, but this is alright (rounding error from floats)
     assert_eq!(a.inverse().unwrap(), i);
}

#[test]
fn mat4_inverse_mult() {
     let a: Matrix4D = [  3.0, -9.0,  7.0,  3.0, 
                          3.0,  8.0,  2.0, -9.0,
                         -4.0,  4.0,  4.0,  1.0,
                         -6.0,  5.0, -1.0,  1.0, ].into();

     let b: Matrix4D = [ 8.0,  2.0, 2.0, 2.0,
                         3.0, -1.0, 7.0, 0.0,
                         7.0,  0.0, 5.0, 4.0,
                         6.0, -2.0, 0.0, 5.0  ].into();

     let c = a * b;

     // Fails, but this is alright (rounding error from floats)
     assert_eq!(a, c * b.inverse().unwrap());
}

#[test]
fn mat4_translation() {
    let transform = Matrix4D::translation(5.0, -3.0, 2.0);
    let point = Tuple4D::point(-3.0, 4.0, 5.0);

    assert_eq!(transform * point, Tuple4D::point(2.0, 1.0, 7.0));
}

#[test]
fn mat4_translation_inverse() {
    let transform = Matrix4D::translation(5.0, -3.0, 2.0).inverse().unwrap();
    let point = Tuple4D::point(-3.0, 4.0, 5.0);

    assert_eq!(transform * point, Tuple4D::point(-8.0, 7.0, 3.0));
}

#[test]
fn mat4_translation_vector() {
    let transform = Matrix4D::translation(5.0, -3.0, 2.0);
    let vector = Tuple4D::vector(-3.0, 4.0, 5.0);

    assert_eq!(transform * vector, vector);
}

#[test]
fn mat4_scaling() {
    let transform = Matrix4D::scaling(2.0, 3.0, 4.0);
    let vector = Tuple4D::vector(-4.0, 6.0, 8.0);

    assert_eq!(transform * vector, Tuple4D::vector(-8.0, 18.0, 32.0));
}

#[test]
fn mat4_scaling_inverse() {
    let transform = Matrix4D::scaling(2.0, 3.0, 4.0).inverse().unwrap();
    let vector = Tuple4D::vector(-4.0, 6.0, 8.0);

    assert_eq!(transform * vector, Tuple4D::vector(-2.0, 2.0, 2.0));
}

#[test]
fn mat4_scaling_reflection() {
    let transform = Matrix4D::scaling(-1.0, 1.0, 1.0);
    let point = Tuple4D::point(2.0, 3.0, 4.0);

    assert_eq!(transform * point, Tuple4D::point(-2.0, 3.0, 4.0));
}

#[test]
fn mat4_rotate_x() {
    let half_quarter = Matrix4D::rotation_x(std::f64::consts::PI / 4.0);
    let full_quarter = Matrix4D::rotation_x(std::f64::consts::PI / 2.0);
    let point = Tuple4D::point(0.0, 1.0, 0.0);

    assert_eq!(full_quarter * point,
        Tuple4D::point(0.0, 0.0, 1.0));
    assert_eq!(half_quarter * point,
        Tuple4D::point(0.0, 2.0f64.sqrt() / 2.0, 2.0f64.sqrt() / 2.0));
}

#[test]
fn mat4_rotate_y() {
    let half_quarter = Matrix4D::rotation_y(std::f64::consts::PI / 4.0);
    let full_quarter = Matrix4D::rotation_y(std::f64::consts::PI / 2.0);
    let point = Tuple4D::point(0.0, 0.0, 1.0);

    assert_eq!(full_quarter * point,
        Tuple4D::point(1.0, 0.0, 0.0));
    assert_eq!(half_quarter * point,
        Tuple4D::point(2.0f64.sqrt() / 2.0, 0.0, 2.0f64.sqrt() / 2.0));
}

#[test]
fn mat4_rotate_z() {
    let half_quarter = Matrix4D::rotation_z(std::f64::consts::PI / 4.0);
    let full_quarter = Matrix4D::rotation_z(std::f64::consts::PI / 2.0);
    let point = Tuple4D::point(0.0, 1.0, 0.0);

    assert_eq!(full_quarter * point,
        Tuple4D::point(-1.0, 0.0, 0.0));
    assert_eq!(half_quarter * point,
        Tuple4D::point(-2.0f64.sqrt() / 2.0, 2.0f64.sqrt() / 2.0, 0.0));
}

#[test]
fn mat4_shear_xy() {
    let transform = Matrix4D::shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let point = Tuple4D::point(2.0, 3.0, 4.0);

    assert_eq!(transform * point, Tuple4D::point(5.0, 3.0, 4.0));
}

#[test]
fn mat4_shear_xz() {
    let transform = Matrix4D::shearing(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    let point = Tuple4D::point(2.0, 3.0, 4.0);

    assert_eq!(transform * point, Tuple4D::point(6.0, 3.0, 4.0));
}

#[test]
fn mat4_shear_yx() {
    let transform = Matrix4D::shearing(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    let point = Tuple4D::point(2.0, 3.0, 4.0);

    assert_eq!(transform * point, Tuple4D::point(2.0, 5.0, 4.0));
}

#[test]
fn mat4_shear_yz() {
    let transform = Matrix4D::shearing(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    let point = Tuple4D::point(2.0, 3.0, 4.0);

    assert_eq!(transform * point, Tuple4D::point(2.0, 7.0, 4.0));
}

#[test]
fn mat4_shear_zx() {
    let transform = Matrix4D::shearing(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    let point = Tuple4D::point(2.0, 3.0, 4.0);

    assert_eq!(transform * point, Tuple4D::point(2.0, 3.0, 6.0));
}

#[test]
fn mat4_shear_zy() {
    let transform = Matrix4D::shearing(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    let point = Tuple4D::point(2.0, 3.0, 4.0);

    assert_eq!(transform * point, Tuple4D::point(2.0, 3.0, 7.0));
}

#[test]
fn chained_transforms() {
    let a = Matrix4D::rotation_x(std::f64::consts::PI / 2.0);
    let b = Matrix4D::scaling(5.0, 5.0, 5.0);
    let c = Matrix4D::translation(10.0, 5.0, 7.0);

    let t = c * b * a;
    let p = Tuple4D::point(1.0, 0.0, 1.0);

    assert_eq!(t * p, Tuple4D::point(15.0, 0.0, 7.0));
}

#[test]
fn default_view() {
    let from = Tuple4D::point(0.0, 0.0, 0.0);
    let to = Tuple4D::point(0.0, 0.0, -1.0);
    let up = Tuple4D::vector(0.0, 1.0, 0.0);

    assert_eq!(Matrix4D::identity(), Matrix4D::view_transform(from, to, up));
}

#[test]
fn positive_z_view() {
    let from = Tuple4D::point(0.0, 0.0, 0.0);
    let to = Tuple4D::point(0.0, 0.0, 1.0);
    let up = Tuple4D::vector(0.0, 1.0, 0.0);

    assert_eq!(Matrix4D::view_transform(from, to, up),
        Matrix4D::scaling(-1.0, 1.0, -1.0));
}

#[test]
fn view_moves_world() {
    let from = Tuple4D::point(0.0, 0.0, 8.0);
    let to = Tuple4D::point(0.0, 0.0, 0.0);
    let up = Tuple4D::vector(0.0, 1.0, 0.0);

    assert_eq!(Matrix4D::view_transform(from, to, up),
        Matrix4D::translation(0.0, 0.0, -8.0));
}

#[test]
fn arbitrary_view() {
    let from = Tuple4D::point(1.0, 3.0, 2.0);
    let to = Tuple4D::point(4.0, -2.0, 8.0);
    let up = Tuple4D::vector(1.0, 1.0, 0.0);

    let a: Matrix4D = [  -0.50709, 0.50709,  0.67612, -2.36643, 
                          0.76772, 0.60609,  0.12122, -2.82843,
                         -0.35857, 0.59761, -0.71714,  0.00000,
                         -0.00000, 0.00000,  0.00000,  1.00000, ].into();

    assert_eq!(Matrix4D::view_transform(from, to, up), a);
}
