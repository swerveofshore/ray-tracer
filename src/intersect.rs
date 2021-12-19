use crate::consts::FEQ_EPSILON;
use crate::tuple::Tuple4D;
use crate::ray::Ray4D;
use crate::geometry::{ ShapeDebug, normal_at };

/// An intersection.
///
/// This structure assumes that some ray produced an intersection. Parameter `t`
/// is analogous to `t` for a ray (the offset from the ray origin).
///
/// The `what` parameter is a reference to an `IntersectableDebug` trait object.
/// A trait object is used because there can be many unique objects (sphere,
/// cube, plane, etc.) which are intersectable.
#[derive(Copy, Clone, Debug)]
pub struct Intersection<'a> {
    pub t: f64,
    pub what: &'a dyn ShapeDebug,
}

/// Implements partial equality on an Intersection.
///
/// Two Intersection structures are equal if the offsets `t` of the
/// intersections are equivalent, and if the underlying *pointers* of the
/// intersections are equivalent.
impl<'a> PartialEq for Intersection<'a> {
    fn eq(&self, other: &Intersection<'a>) -> bool {
        self.t == other.t && std::ptr::eq(self.what, other.what)
    }
}

/// A collection of intersections.
///
/// Mostly a wrapper for a vector of `Intersection` objects. See the
/// `Intersection` documentation for more information.
#[derive(Clone, Debug, Default)]
pub struct Intersections<'a> {
    pub intersections: Vec<Intersection<'a>>,
}

impl<'a> Intersections<'a> {
    /// Creates a new list of intersections.
    pub fn new() -> Intersections<'a> {
        Intersections { intersections: Vec::new() }
    }

    /// Checks if any Intersectable object has been hit.
    ///
    /// If no hit is registered, this function returns `None`.
    ///
    /// Effectively, a hit occurs if at least one `Intersection` in this
    /// `Intersections` is finite and is greater than or equal to 0.
    ///
    /// As a note, this function sorts the `intersections` field on every call.
    /// This is because the `Intersection` with the lowest `t` is chosen. A more
    /// optimal implementation is likely possible.
    pub fn hit<'b>(&'b mut self) -> Option<Intersection<'a>> {
        self.intersections.retain(|i| i.t.is_finite());
        self.sort();

        for i in self.intersections.iter() {
            if i.t >= 0.0 {
                return Some(*i);
            }
        }

        None
    }

    /// Sorts the intersections by `t`, ignoring `f64` semantics.
    pub fn sort(&mut self) {
        self.intersections.sort_by(|a, b|
            a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal)
        );
    }
}

/// Intersects a ray with a `Shape`.
///
/// Each shape implements `local_intersect`, which calculates intersections of a
/// shape in *normal* space. This function uses the `transform` associated with
/// each `Shape`, and converts the ray to object space before calculating
/// intersections.
///
/// Intersections are technically calculated in local space; these local-space
/// intersections are returned in an `Intersections` record.
pub fn intersect(s: & dyn ShapeDebug, r: Ray4D) -> Intersections {
    let inverse_transform = s.transform().inverse().expect(
        "Transformation matrix on shape should be invertible."
    );

    let transformed_ray = r.transform(inverse_transform);
    s.local_intersect(transformed_ray)
}

/// A record for computations associated with an `Intersection`.
///
/// Mostly a superset of an `Intersection`.
#[derive(Clone, Debug)]
pub struct IntersectionComputation<'a> {
    /// The "time" of the ray intersection.
    pub t: f64,

    /// The object being intersected.
    pub obj: &'a dyn ShapeDebug,

    /// The point where the intersection occurs.
    pub point: Tuple4D,

    /// A point slightly above the intersected surface. Used to prevent an
    /// object from shadowing itself (this causes "acne").
    pub over_point: Tuple4D,

    /// A point slightly below the intersected surface. Used to prevent an
    /// object from refracting itself on its surface.
    pub under_point: Tuple4D,

    /// The eye vector for the intersection.
    pub eyev: Tuple4D,

    /// The normal vector of the object being intersected.
    pub normalv: Tuple4D,

    /// The intersection ray, reflected across the normal.
    pub reflectv: Tuple4D,

    /// Whether the intersection occurs within the object or not.
    pub inside: bool,

    /// The refractive index of the material being exited.
    pub n1: f64,

    /// The refractive index of the material being entered.
    pub n2: f64,
}

impl<'a> IntersectionComputation<'a> {
    /// Creates a new intersection computation, given a ray and intersection.
    ///
    /// The `is` parameter is a collection of intersections. If provided,
    /// refraction indices will be calculated.
    pub fn new(r: &Ray4D, hit: &'a Intersection, is: Option<&Intersections<'a>>)
        -> IntersectionComputation<'a> {
        let t = hit.t;
        let obj = hit.what;
        let point = r.position(t);
        let eyev = -r.direction;
        let mut normalv = normal_at(obj, point);

        let inside = if normalv.dot(&eyev) < 0.0 {
            normalv = -normalv;
            true
        } else {
            false
        };

        let over_point = point + normalv * FEQ_EPSILON;
        let under_point = point - normalv * FEQ_EPSILON;

        let reflectv = r.direction.reflect(&normalv);
        let (n1, n2) = if let Some(xs) = is {
            Self::refraction_indices(hit, xs)
        } else {
            (1.0, 1.0)
        };

        IntersectionComputation {
            t, obj,
            point, over_point, under_point,
            eyev, normalv, reflectv,
            inside,
            n1, n2,
        }
    }

    fn refraction_indices(hit: &'a Intersection, is: &Intersections<'a>)
        -> (f64, f64) {
        // The exiting and entering refractive indices, respectively.
        let mut n1 = 1.0;
        let mut n2 = 1.0;

        // Contains all objects which have been encountered, but not yet exited
        // by the refracting ray.
        let mut containers: Vec<&'a dyn ShapeDebug> = Vec::new();

        for i in is.intersections.iter() {
            if i == hit {
                if containers.is_empty() {
                    n1 = 1.0; 
                } else {
                    n1 = containers.last().unwrap().material().refractive_index;
                }
            }

            // If object `i.what` is in `containers`, remove it.
            if let Some(j) 
                = containers.iter().position(|&x| std::ptr::eq(x, i.what)) {
                containers.remove(j);
            }
            // Otherwise, add the object to `containers`.
            else {
                containers.push(i.what);
            }

            if i == hit {
                if containers.is_empty() {
                    n2 = 1.0;
                } else {
                    n2 = containers.last().unwrap().material().refractive_index;
                }

                // TODO: refactor this to be more Rust-idiomatic; the hanging
                // break on this second if feels awkward
                break;
            }
        }

        (n1, n2)
    }

    /// Calculates the reflectance of a hit.
    ///
    /// The reflectance is a number between 0 and 1, representing what fraction
    /// of the light is reflected for the hit.
    ///
    /// TODO: Explain why the below calculations work.
    pub fn schlick(&self) -> f64 {
        let mut cos = self.eyev.dot(&self.normalv);

        // Total internal reflection can only occur if n1 > n2.
        if self.n1 > self.n2 {
            let n = self.n1 / self.n2;
            let sin2_t = n.powi(2) * (1.0 - cos.powi(2));

            // If the squared sin of theta t is greater than 1, return full
            // reflectance.
            if sin2_t > 1.0 {
                return 1.0
            }

            // Otherwise, continue calculations.
            cos = (1.0 - sin2_t).sqrt();
        }

        let r0 = ((self.n1 - self.n2) / (self.n1 + self.n2)).powi(2);
        r0 + (1.0 - r0) * (1.0 - cos).powi(5)
    }
}
