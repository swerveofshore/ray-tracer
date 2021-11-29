pub mod tuple;
pub mod matrix;
pub mod ray;
pub mod light;

pub mod geometry;
pub mod world;
pub mod camera;

pub mod color;
pub mod canvas;

pub mod extra;

const FEQ_EPSILON: f64 = 0.0001;
pub fn feq(left: f64, right: f64) -> bool {
    (left - right).abs() < FEQ_EPSILON
}
