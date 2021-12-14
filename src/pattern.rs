use crate::feq;
use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::color::Color;
use crate::geometry::ShapeDebug;

#[derive(Copy, Clone, Debug, PartialEq)]
enum PatternType {
    /// Does nothing. Pixels are black at every point. Mostly for testing.
    Null,

    /// Returns a color equal to the supplied point.
    Point,

    /// Returns the same color at all points. Mostly for testing.
    Identity(Color),

    /// Creates an stripe pattern, alternating colors across the X axis.
    Stripe(Color, Color),
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Pattern {
    ty: PatternType,
    pub transform: Matrix4D,
}

impl Pattern {
    pub fn null() -> Pattern {
        Pattern {
            ty: PatternType::Null,
            transform: Matrix4D::identity(),
        }
    }

    pub fn point() -> Pattern {
        Pattern {
            ty: PatternType::Point,
            transform: Matrix4D::identity(),
        }
    }

    pub fn identity(c: Color) -> Pattern {
        Pattern {
            ty: PatternType::Identity(c),
            transform: Matrix4D::identity(),
        }
    }

    pub fn stripe(primary: Color, secondary: Color)
        -> Pattern {
        Pattern {
            ty: PatternType::Stripe(primary, secondary),
            transform: Matrix4D::identity(),
        }
    }

    pub fn pattern_at(&self, p: Tuple4D) -> Color {
        match self.ty {
            PatternType::Null => Color::black(),
            PatternType::Point => Color::rgb(p.x, p.y, p.z),
            PatternType::Identity(c) => c,
            PatternType::Stripe(_, _) => self.stripe_at(p),
        }
    }

    fn stripe_at(&self, p: Tuple4D) -> Color {
        // The pattern type here can't be anything *but* a Stripe.
        let (primary, secondary) =
            if let PatternType::Stripe(pr, se) = self.ty {
                (pr, se)
            } else {
                unreachable!();
            };

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
    pub fn pattern_at_object(&self, object: & dyn ShapeDebug, p: Tuple4D)
        -> Color {
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

#[test]
fn point_pattern_with_object_transformation() {
    use crate::geometry::Sphere;

    let mut s = Sphere::unit();
    s.transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let p = Pattern::point();
    assert_eq!(p.pattern_at_object(&s, Tuple4D::point(2.0, 3.0, 4.0)),
        Color::rgb(1.0, 1.5, 2.0));
}

#[test]
fn point_pattern_with_pattern_transformation() {
    use crate::geometry::Sphere;

    let s = Sphere::unit();

    let mut p = Pattern::point();
    p.transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    assert_eq!(p.pattern_at_object(&s, Tuple4D::point(2.0, 3.0, 4.0)),
        Color::rgb(1.0, 1.5, 2.0));
}

#[test]
fn point_pattern_with_object_and_pattern_transformation() {
    use crate::geometry::Sphere;

    let mut s = Sphere::unit();
    s.transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let mut p = Pattern::point();
    p.transform = Matrix4D::translation(0.5, 1.0, 1.5);

    assert_eq!(p.pattern_at_object(&s, Tuple4D::point(2.5, 3.0, 3.5)),
        Color::rgb(0.75, 0.5, 0.25));
}
