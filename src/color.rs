use std::ops::{ Add, Sub, Mul };

use crate::feq;

#[derive(Copy, Clone, Debug, Default, PartialOrd)]
pub struct Color {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl PartialEq for Color {
    fn eq(&self, other: &Color) -> bool {
        feq(self.r, other.r) &&
            feq(self.g, other.g) &&
            feq(self.b, other.b)
    }
}

impl Color {
    pub fn rgb(r: f64, g: f64, b: f64) -> Color {
        Color {
            r, // : r.clamp(0.0, 1.0),
            g, // : g.clamp(0.0, 1.0),
            b, // : b.clamp(0.0, 1.0),
        }
    }

    pub fn black() -> Color {
        Color {
            r: 0.0,
            g: 0.0,
            b: 0.0
        }
    }

    pub fn white() -> Color {
        Color {
            r: 1.0,
            g: 1.0,
            b: 1.0
        }
    }

    pub fn red() -> Color {
        Color {
            r: 1.0,
            g: 0.0,
            b: 0.0
        }
    }

    pub fn green() -> Color {
        Color {
            r: 0.0,
            g: 1.0,
            b: 0.0
        }
    }

    pub fn blue() -> Color {
        Color {
            r: 0.0,
            g: 0.0,
            b: 1.0
        }
    }

    pub fn hadamard(c1: &Color, c2: &Color) -> Color {
        let r = c1.r * c2.r;
        let g = c1.g * c2.g;
        let b = c1.b * c2.b;

        Color { r, g, b }
    }

    pub fn average(c1: &Color, c2: &Color) -> Color {
        let r = (c1.r + c2.r) / 2.0;
        let g = (c1.g + c2.g) / 2.0;
        let b = (c1.b + c2.b) / 2.0;

        Color { r, g, b }
    }
}

impl Add<Color> for Color {
    type Output = Color;

    fn add(self, other: Color) -> Self::Output {
        Color {
            r: self.r + other.r,
            g: self.g + other.g,
            b: self.b + other.b,
        }
    }
}

impl Sub<Color> for Color {
    type Output = Color;

    fn sub(self, other: Color) -> Self::Output {
        Color {
            r: self.r - other.r,
            g: self.g - other.g,
            b: self.b - other.b,
        }
    }
}

impl Mul<f64> for Color {
    type Output = Color;

    fn mul(self, other: f64) -> Self::Output {
        Color {
            r: self.r * other,
            g: self.g * other,
            b: self.b * other,
        }
    }
}

impl Mul<Color> for f64 {
    type Output = Color;

    fn mul(self, other: Color) -> Self::Output {
        Color {
            r: self * other.r,
            g: self * other.g,
            b: self * other.b,
        }
    }
}

impl Mul<Color> for Color {
    type Output = Color;

    fn mul(self, other: Color) -> Self::Output {
        Color::hadamard(&self, &other)
    }
}

#[test]
fn add_colors() {
    let c1 = Color::rgb(0.9, 0.6, 0.75);
    let c2 = Color::rgb(0.7, 0.1, 0.25);
    let c3 = Color { r: 1.6, g: 0.7, b: 1.0 };

    assert_eq!(c1 + c2, c3);
}

#[test]
fn sub_colors() {
    let c1 = Color::rgb(0.9, 0.6, 0.75);
    let c2 = Color::rgb(0.7, 0.1, 0.25);
    let c3 = Color { r: 0.2, g: 0.5, b: 0.5 };

    assert_eq!(c1 - c2, c3);
}

#[test]
fn mul_color() {
    let c1 = Color::rgb(0.2, 0.3, 0.4);
    let c2 = Color { r: 0.4, g: 0.6, b: 0.8 };

    assert_eq!(c1 * 2.0, c2);
}
