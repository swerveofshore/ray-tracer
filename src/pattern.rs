use crate::feq;
use crate::tuple::Tuple4D;
use crate::color::Color;

pub trait Pattern {
    fn pattern_at(&self, p: Tuple4D) -> Color;
}

/// An alternating stripe pattern applied across the X axis.
///
/// Effectively, for a point `(x, y, z)`, if `floor(x) % 2 == 0`, the `primary'
/// color is applied to that point; otherwise, the `secondary` color is used.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct StripePattern {
    primary: Color,
    secondary: Color,
}

impl StripePattern {
    pub fn new(primary: Color, secondary: Color) -> StripePattern {
        StripePattern { primary, secondary }
    }

    pub fn stripe_at(&self, p: Tuple4D) -> Color {
        if feq(p.x.floor().rem_euclid(2.0), 0.0) {
            self.primary
        } else {
            self.secondary
        }
    }
}

impl Pattern for StripePattern {
    fn pattern_at(&self, p: Tuple4D) -> Color {
        self.stripe_at(p)
    }
}

#[test]
fn stripe_pattern_is_constant_along_y() {
    let pattern = StripePattern::new(Color::white(), Color::black());

    assert_eq!(pattern.stripe_at(Tuple4D::point(0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.stripe_at(Tuple4D::point(0.0, 1.0, 0.0)),
        Color::white());
    assert_eq!(pattern.stripe_at(Tuple4D::point(0.0, 2.0, 0.0)),
        Color::white());
}

#[test]
fn stripe_pattern_is_constant_along_z() {
    let pattern = StripePattern::new(Color::white(), Color::black());

    assert_eq!(pattern.stripe_at(Tuple4D::point(0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.stripe_at(Tuple4D::point(0.0, 0.0, 1.0)),
        Color::white());
    assert_eq!(pattern.stripe_at(Tuple4D::point(0.0, 0.0, 2.0)),
        Color::white());
}

#[test]
fn stripe_pattern_alternates_along_x() {
    let pattern = StripePattern::new(Color::white(), Color::black());

    assert_eq!(pattern.stripe_at(Tuple4D::point( 0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.stripe_at(Tuple4D::point( 0.9, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.stripe_at(Tuple4D::point( 1.0, 0.0, 0.0)),
        Color::black());
    assert_eq!(pattern.stripe_at(Tuple4D::point(-0.1, 0.0, 0.0)),
        Color::black());
    assert_eq!(pattern.stripe_at(Tuple4D::point(-1.0, 0.0, 0.0)),
        Color::black());
    assert_eq!(pattern.stripe_at(Tuple4D::point(-1.1, 0.0, 0.0)),
        Color::white());
}
