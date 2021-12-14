use crate::feq;
use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::color::Color;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum PatternType {
    Stripe(Color, Color),
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Pattern {
    ty: PatternType,
    transform: Matrix4D,
}

impl Pattern {
    pub fn stripe(transform: Matrix4D, primary: Color, secondary: Color)
        -> Pattern {
        Pattern { ty: PatternType::Stripe(primary, secondary), transform }
    }

    pub fn pattern_at(&self, p: Tuple4D) -> Color {
        match self.ty {
            PatternType::Stripe(_, _) => self.stripe_at(p),
        }
    }

    fn stripe_at(&self, p: Tuple4D) -> Color {
        // The pattern type here can't be anything *but* a Stripe.
        let PatternType::Stripe(primary, secondary) = self.ty;

        if feq(p.x.floor().rem_euclid(2.0), 0.0) {
            primary
        } else {
            secondary
        }
    }
}

#[test]
fn stripe_pattern_is_constant_along_y() {
    let pattern = Pattern::stripe(Matrix4D::identity(),
        Color::white(), Color::black());

    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 1.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 2.0, 0.0)),
        Color::white());
}

#[test]
fn stripe_pattern_is_constant_along_z() {
    let pattern = Pattern::stripe(Matrix4D::identity(),
        Color::white(), Color::black());

    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 1.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 2.0)),
        Color::white());
}

#[test]
fn stripe_pattern_alternates_along_x() {
    let pattern = Pattern::stripe(Matrix4D::identity(),
        Color::white(), Color::black());

    assert_eq!(pattern.pattern_at(Tuple4D::point( 0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point( 0.9, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point( 1.0, 0.0, 0.0)),
        Color::black());
    assert_eq!(pattern.pattern_at(Tuple4D::point(-0.1, 0.0, 0.0)),
        Color::black());
    assert_eq!(pattern.pattern_at(Tuple4D::point(-1.0, 0.0, 0.0)),
        Color::black());
    assert_eq!(pattern.pattern_at(Tuple4D::point(-1.1, 0.0, 0.0)),
        Color::white());
}
