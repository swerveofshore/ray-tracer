use crate::feq;
use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::color::Color;
use crate::geometry::ShapeDebug;

#[derive(Copy, Clone, Debug, PartialEq)]
enum PatternType {
    Stripe(Color, Color),
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Pattern {
    ty: PatternType,
    transform: Matrix4D,
}

impl Pattern {
    pub fn stripe(primary: Color, secondary: Color)
        -> Pattern {
        Pattern {
            ty: PatternType::Stripe(primary, secondary),
            transform: Matrix4D::identity(),
        }
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

    /// Applies a Pattern to a Shape.
    ///
    /// Assumes that point `p` is in world space. Before a pattern is applied to
    /// point `p`, `p` is transformed to object space (local to the `Shape`),
    /// afterwards transformed to pattern space (local to the `Pattern` on the
    /// `Shape`).
    fn pattern_at_object(&self, object: & dyn ShapeDebug, p: Tuple4D) -> Color {
        let oti = object.transform().inverse().expect(
            "Object transform should be invertible."
        );

        let pti = self.transform.inverse().expect(
            "Pattern transform should be invertible."
        );

        let object_point = oti * p;
        let pattern_point = pti * object_point;

        self.pattern_at(pattern_point)
    }
}

#[test]
fn stripe_pattern_is_constant_along_y() {
    let pattern = Pattern::stripe(Color::white(), Color::black());

    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 1.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 2.0, 0.0)),
        Color::white());
}

#[test]
fn stripe_pattern_is_constant_along_z() {
    let pattern = Pattern::stripe(Color::white(), Color::black());

    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 0.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 1.0)),
        Color::white());
    assert_eq!(pattern.pattern_at(Tuple4D::point(0.0, 0.0, 2.0)),
        Color::white());
}

#[test]
fn stripe_pattern_alternates_along_x() {
    let pattern = Pattern::stripe(Color::white(), Color::black());

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

#[test]
fn stripes_on_object_transformation() {
    use crate::geometry::Sphere;

    let mut object = Sphere::unit();
    object.transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let pattern = Pattern::stripe(Color::white(), Color::black());

    let p = Tuple4D::point(1.5, 0.0, 0.0);
    assert_eq!(pattern.pattern_at_object(&object, p), Color::white());
}

#[test]
fn stripes_on_pattern_transformation() {
    use crate::geometry::Sphere;

    let object = Sphere::unit();

    let mut pattern = Pattern::stripe(Color::white(), Color::black());
    pattern.transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let p = Tuple4D::point(1.5, 0.0, 0.0);
    assert_eq!(pattern.pattern_at_object(&object, p), Color::white());
}

#[test]
fn stripes_on_object_and_pattern_transformation() {
    use crate::geometry::Sphere;

    let mut object = Sphere::unit();
    object.transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let mut pattern = Pattern::stripe(Color::white(), Color::black());
    pattern.transform = Matrix4D::translation(0.5, 0.0, 0.0);

    let p = Tuple4D::point(2.5, 0.0, 0.0);
    assert_eq!(pattern.pattern_at_object(&object, p), Color::white());
}
