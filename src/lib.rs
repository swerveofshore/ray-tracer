pub mod tuple;
pub mod matrix;

pub mod ray;
pub mod geometry;
pub mod intersect;

pub mod light;
pub mod color;
pub mod pattern;

pub mod world;
pub mod camera;
pub mod canvas;

pub mod extra;
pub mod consts;

pub fn feq(left: f64, right: f64) -> bool {
    (left - right).abs() < consts::FEQ_EPSILON
}
