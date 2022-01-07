use std::ops::{ Add, Sub, Mul };

use crate::feq;

/// A color.
///
/// Represented conventionally with red-green-blue (RGB) values. Each value
/// ranges from 0.0 to 1.0 inclusive.
///
/// # Examples
///
/// Construct the color red:
///
/// ```
/// # #![allow(unused)]
/// # use ray_tracer_challenge::color::Color;
/// let red = Color::red();
/// assert_eq!(red, Color::rgb(1.0, 0.0, 0.0));
/// ```
///
/// Blend two colors:
///
/// ```
/// # #![allow(unused)]
/// # use ray_tracer_challenge::color::Color;
/// let green = Color::green();
/// let blue = Color::blue();
/// let blend = Color::average(&green, &blue);
/// assert_eq!(blend, Color::rgb(0.0, 0.5, 0.5));
/// ```
#[derive(Copy, Clone, Debug, Default, PartialOrd)]
pub struct Color {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

/// Partial equality on two colors.
///
/// Similar to the `PartialEq` implementation on `Tuple4D`, `Color`s are
/// compared component-wise, accounting for possible floating point error in
/// comparisons.
impl PartialEq for Color {
    fn eq(&self, other: &Color) -> bool {
        feq(self.r, other.r) &&
            feq(self.g, other.g) &&
            feq(self.b, other.b)
    }
}

/// Conversion from a vector to a `Color`.
///
/// Takes the first `n` elements of a vector, and assigns them to the `r`, `g`
/// and `b` fields of the `Color`, in that order. If there aren't enough
/// elements in the vector (e.g. `n == 2`), fields are assigned defaults in
/// place.
///
/// For instance, for a vector `v` of length `2`, `r == v[0]` and `g == v[1]`.
impl From<&Vec<f64>> for Color {
    fn from(v: &Vec<f64>) -> Color {
        match v.len() {
            0 => Default::default(),
            1 => Color { r: v[0], ..Default::default() },
            2 => Color { r: v[0], g: v[1], ..Default::default() },
            _ => Color { r: v[0], g: v[1], b: v[2] }
        }
    }
}

impl Color {
    /// Creates a color with red, green and blue values.
    pub fn rgb(r: f64, g: f64, b: f64) -> Color {
        Color { r, g, b }
    }

    /// The color black.
    pub fn black() -> Color {
        Color {
            r: 0.0,
            g: 0.0,
            b: 0.0
        }
    }

    /// The color white.
    pub fn white() -> Color {
        Color {
            r: 1.0,
            g: 1.0,
            b: 1.0
        }
    }

    /// The color red.
    pub fn red() -> Color {
        Color {
            r: 1.0,
            g: 0.0,
            b: 0.0
        }
    }

    /// The color green.
    pub fn green() -> Color {
        Color {
            r: 0.0,
            g: 1.0,
            b: 0.0
        }
    }

    /// The color blue.
    pub fn blue() -> Color {
        Color {
            r: 0.0,
            g: 0.0,
            b: 1.0
        }
    }

    /// Computes the Hadamard product of two colors.
    ///
    /// This is provided as an associated function of `Color` to prevent
    /// possible confusion with the `*` operator.
    ///
    /// The hadamard product multiplies each component of the two colors, and
    /// yields a new color containing those products. (In other words, this is
    /// a dot product which doesn't sum the component products).
    ///
    /// # Examples
    ///
    /// Computing the Hadamard product between yellow and purple:
    ///
    /// ```
    /// # use ray_tracer_challenge::color::Color;
    /// let yellow = Color::rgb(1.0, 1.0, 0.0);
    /// let purple = Color::rgb(1.0, 0.0, 1.0);
    /// let product = Color::hadamard(&yellow, &purple);
    /// assert_eq!(product, Color::red());
    /// ```
    pub fn hadamard(c1: &Color, c2: &Color) -> Color {
        let r = c1.r * c2.r;
        let g = c1.g * c2.g;
        let b = c1.b * c2.b;

        Color { r, g, b }
    }

    /// Averages two colors.
    ///
    /// This is provided as an associated function of `Color` to prevent
    /// possible confusion with the `*` operator.
    ///
    /// Each component of a color is averaged with the other color's respective
    /// component. For example, for colors `c1` and `c2`, the average color's
    /// red component is equal to `(c1.r + c2.r) / 2.0`.
    ///
    /// # Examples
    ///
    /// Computing the average between cyan and purple:
    ///
    /// ```
    /// # use ray_tracer_challenge::color::Color;
    /// let cyan = Color::rgb(0.0, 1.0, 1.0);
    /// let purple = Color::rgb(1.0, 0.0, 1.0);
    /// let avg = Color::average(&cyan, &purple);
    /// assert_eq!(avg, Color::rgb(0.5, 0.5, 1.0));
    /// ```
    pub fn average(c1: &Color, c2: &Color) -> Color {
        let r = (c1.r + c2.r) / 2.0;
        let g = (c1.g + c2.g) / 2.0;
        let b = (c1.b + c2.b) / 2.0;

        Color { r, g, b }
    }
}

/// Adds two colors together.
///
/// Components are added together individually.
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

/// Subtracts one color from another.
///
/// Components are subtracted from one another individually.
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

/// Multiplies a color by a scalar.
///
/// Each component is multiplied by the scalar.
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

/// Multiplies a scalar by a color.
///
/// Returns a color with each component multiplied by the scalar.
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

/// Multiplies a color by a color.
///
/// For colors `c1` and `c2`, `c1 * c2` is shorthand for
/// `Color::hadamard(&c1, &c2)`. Consider using the `Color::hadamard` call
/// directly, as the `*` operator for `Color`s can be somewhat ambiguous.
///
/// # Examples
///
/// ```
/// # use ray_tracer_challenge::color::Color;
/// let c1 = Color::red();
/// let c2 = Color::blue();
/// assert_eq!(c1 * c2, Color::hadamard(&c1, &c2));
/// ```
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
fn subtract_colors() {
    let c1 = Color::rgb(0.9, 0.6, 0.75);
    let c2 = Color::rgb(0.7, 0.1, 0.25);
    let c3 = Color { r: 0.2, g: 0.5, b: 0.5 };

    assert_eq!(c1 - c2, c3);
}

#[test]
fn multiply_colors() {
    let c1 = Color::rgb(0.2, 0.3, 0.4);
    let c2 = Color { r: 0.4, g: 0.6, b: 0.8 };

    assert_eq!(c1 * 2.0, c2);
}
