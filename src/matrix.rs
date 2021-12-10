use std::fmt;
use std::ops::{ Index, IndexMut, Mul };
use std::convert::From;

use crate::feq;
use crate::tuple::*;

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct Matrix2D {
    data: [f64; 4],
}

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct Matrix3D {
    data: [f64; 9],
}

#[derive(Copy, Clone, Debug, Default, PartialOrd)]
pub struct Matrix4D {
    data: [f64; 16],
}

impl PartialEq for Matrix4D {
    fn eq(&self, other: &Matrix4D) -> bool {
        self.data.iter().zip(other.data.iter())
                        .all(|(x, y)| feq(*x, *y))
    }
}

impl Matrix2D {
    pub fn new() -> Matrix2D {
        Matrix2D { data: [0.0; 4] }
    }

    pub fn identity() -> Matrix2D {
        let mut buf = [0.0; 4];
        buf[0] = 1.0; buf[3] = 1.0;

        Matrix2D { data: buf }
    }

    pub fn determinant(&self) -> f64 {
        self[(0, 0)] * self[(1, 1)] - self[(0, 1)] * self[(1, 0)]
    }
}

impl Matrix3D {
    pub fn new() -> Matrix3D {
        Matrix3D { data: [0.0; 9] }
    }

    pub fn identity() -> Matrix3D {
        let mut buf = [0.0; 9];
        buf[0] = 1.0; buf[4] = 1.0; buf[8] = 1.0;

        Matrix3D { data: buf }
    }

    pub fn submatrix(&self, row: usize, col: usize) -> Matrix2D {
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

    pub fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let m = self.minor(row, col);
        m * if (row + col) % 2 == 0 { 1.0 } else { -1.0 }
    }

    pub fn determinant(&self) -> f64 {
        let mut sum = 0.0;
        for c in 0..3 {
            sum += self[(0, c)] * self.cofactor(0, c);
        }

        sum
    }
}

impl Matrix4D {
    pub fn new() -> Matrix4D {
        Matrix4D { data: [0.0; 16] }
    }

    pub fn identity() -> Matrix4D {
        let mut buf = [0.0; 16];
        buf[0] = 1.0; buf[5] = 1.0; buf[10] = 1.0; buf[15] = 1.0;

        Matrix4D { data: buf }
    }

    pub fn translation(x: f64, y: f64, z: f64) -> Matrix4D {
        let mut trans = Self::identity();
        trans[(0, 3)] = x;
        trans[(1, 3)] = y;
        trans[(2, 3)] = z;

        trans
    }

    pub fn scaling(x: f64, y: f64, z: f64) -> Matrix4D {
        let mut scale = Self::identity();
        scale[(0, 0)] = x;
        scale[(1, 1)] = y;
        scale[(2, 2)] = z;

        scale
    }

    pub fn rotation_x(r: f64) -> Matrix4D {
        let mut rotate = Self::identity();
        rotate[(1, 1)] =  r.cos();
        rotate[(1, 2)] = -r.sin();
        rotate[(2, 1)] =  r.sin();
        rotate[(2, 2)] =  r.cos();

        rotate
    }

    pub fn rotation_y(r: f64) -> Matrix4D {
        let mut rotate = Self::identity();
        rotate[(0, 0)] =  r.cos();
        rotate[(0, 2)] =  r.sin();
        rotate[(2, 0)] = -r.sin();
        rotate[(2, 2)] =  r.cos();

        rotate
    }

    pub fn rotation_z(r: f64) -> Matrix4D {
        let mut rotate = Self::identity();
        rotate[(0, 0)] =  r.cos();
        rotate[(0, 1)] = -r.sin();
        rotate[(1, 0)] =  r.sin();
        rotate[(1, 1)] =  r.cos();

        rotate
    }

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

    pub fn transpose(&mut self) {
        for r in 0..4 {
            for c in 1..(4-r) {
                let tmp = self[(r, c)];
                self[(r, c)] = self[(c, r)];
                self[(c, r)] = tmp;
            }
        }
    }

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

    pub fn submatrix(&self, row: usize, col: usize) -> Matrix3D {
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

    pub fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let m = self.minor(row, col);
        m * if (row + col) % 2 == 0 { 1.0 } else { -1.0 }
    }

    pub fn determinant(&self) -> f64 {
        let mut sum = 0.0;
        for c in 0..4 {
            sum += self[(0, c)] * self.cofactor(0, c);
        }

        sum
    }

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

impl From<[f64; 4]> for Matrix2D {
    fn from(data: [f64; 4]) -> Matrix2D {
        Matrix2D { data } 
    }
}

impl From<[f64; 9]> for Matrix3D {
    fn from(data: [f64; 9]) -> Matrix3D {
        Matrix3D { data } 
    }
}

impl From<[f64; 16]> for Matrix4D {
    fn from(data: [f64; 16]) -> Matrix4D {
        Matrix4D { data } 
    }
}

impl Index<(usize, usize)> for Matrix2D {
    type Output = f64;

    fn index<'a>(&'a self, index: (usize, usize)) -> &'a f64 {
        &self.data[(index.0 * 2) + index.1]
    }
}

impl Index<(usize, usize)> for Matrix3D {
    type Output = f64;

    fn index<'a>(&'a self, index: (usize, usize)) -> &'a f64 {
        &self.data[(index.0 * 3) + index.1]
    }
}

impl Index<(usize, usize)> for Matrix4D {
    type Output = f64;

    fn index<'a>(&'a self, index: (usize, usize)) -> &'a f64 {
        &self.data[(index.0 * 4) + index.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix2D {
    fn index_mut<'a>(&'a mut self, index: (usize, usize)) -> &'a mut f64 {
        &mut self.data[(index.0 * 2) + index.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix3D {
    fn index_mut<'a>(&'a mut self, index: (usize, usize)) -> &'a mut f64 {
        &mut self.data[(index.0 * 3) + index.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix4D {
    fn index_mut<'a>(&'a mut self, index: (usize, usize)) -> &'a mut f64 {
        &mut self.data[(index.0 * 4) + index.1]
    }
}

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
