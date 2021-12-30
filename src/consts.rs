// Runtime parameters
pub const NUM_THREADS: usize = 3;
pub const CANVAS_WIDTH: usize = 125;
pub const CANVAS_HEIGHT: usize = 100;
pub const OUT_FILE: &'static str = "./out.ppm";

// Floating point comparisons
pub const FEQ_EPSILON: f64 = 0.0001;

// Maximum recursion depths
pub const REFLECTION_RECURSION_DEPTH: usize = 5;
pub const REFRACTION_RECURSION_DEPTH: usize = 5;

// Common refraction indices
pub const VACUUM_RI: f64 = 1.0;
pub const AIR_RI: f64 = 1.00029;
pub const WATER_RI: f64 = 1.333;
pub const GLASS_RI: f64 = 1.52;
pub const DIAMOND_RI: f64 = 2.417;
