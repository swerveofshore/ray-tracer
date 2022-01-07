use crate::consts::FEQ_EPSILON;
use crate::tuple::Tuple4D;
use crate::ray::Ray4D;
use crate::shape::{ Shape, ShapeType, normal_at };

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
    pub what: &'a Shape,

    pub uv: Option<(f64, f64)>,
}

impl<'a> Intersection<'a> {
    /// Create a new intersection at `t` with shape `what`.
    pub fn new(t: f64, what: &'a Shape) -> Intersection<'a> {
        Intersection { t, what, uv: None }
    }

    /// Cratea  new intersection at `t` with shape `what`, specifying UV.
    pub fn new_uv(t: f64, what: &'a Shape, u: f64, v: f64)
        -> Intersection<'a> {
        Intersection { t, what, uv: Some((u, v)) }
    }
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

    /// Aggregates several Intersections together.
    ///
    /// Sorts the combined intersections.
    pub fn aggregate<'b>(mut multi_is: Vec<Intersections<'b>>)
        -> Intersections<'b> {
        let mut combined_is = Intersections::new();
        for is in multi_is.iter_mut() {
            combined_is.intersections.append(&mut is.intersections);
        }

        combined_is.sort();
        combined_is
    }
}

/// A record for computations associated with an `Intersection`.
///
/// Mostly a superset of an `Intersection`.
#[derive(Clone, Debug)]
pub struct IntersectionComputation<'a> {
    /// The "time" of the ray intersection.
    pub t: f64,

    /// The object being intersected.
    pub obj: &'a Shape,

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
        let mut normalv = normal_at(obj, point, hit);

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

    /// Find the refraction indices for a hit in a list of intersections.
    ///
    /// This is mostly important for objects within objects; as a ray passes
    /// from one medium to another, different refraction indices become
    /// relevant.
    fn refraction_indices(hit: &'a Intersection, is: &Intersections<'a>)
        -> (f64, f64) {
        // The exiting and entering refractive indices, respectively.
        let mut n1 = 1.0;
        let mut n2 = 1.0;

        // Contains all objects which have been encountered, but not yet exited
        // by the refracting ray.
        let mut containers: Vec<&'a Shape> = Vec::new();

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

/// Checks whether an intersection is allowed for a CSG operation.
///
/// This function assumes that `op.borrow().ty` is either a `Union`,
/// `Intersection` or `Difference`.
///
/// For nomenclature, this function assumes the existence of two shapes, a
/// "left shape" and a "right shape" (these map to the left and right operands
/// of argument `op`). If `which_hit` is true, the intersection hit the left
/// shape, otherwise the intersection hit the right shape.
///
/// Variable `inside_left` is true if there's an intersection within the left
/// shape, and variable `inside_right` if there's an intersection within the
/// right shape.
///
/// This function returns whether an intersection is allowed for the CSG node.
pub fn intersection_allowed(op: &Shape, which_hit: bool, inside_left: bool,
    inside_right: bool) -> bool {
    match op.ty {
        ShapeType::Union(_, _) =>
            (which_hit && !inside_right) || (!which_hit && !inside_left),
        ShapeType::Intersection(_, _) =>
            (which_hit && inside_right) || (!which_hit && inside_left),
        ShapeType::Difference(_, _) =>
            (which_hit && !inside_right) || (!which_hit && inside_left),
        // By default, disallow intersection
        _ => false,
    }
}

/// Filter intersections based on CSG operations.
///
/// While all shapes in CSG operations are still intersected, technically, we
/// only want certain shape intersections to appear in the resultant render.
///
/// This function removes intersections filtered by the CSG operation. This is
/// most easily visualized with something like a `Difference` operation, which
/// explicitly subtracts one shape from another.
pub fn filter_intersections<'a>(op: &Shape, is: &Intersections<'a>)
    -> Intersections<'a> {
    // Start outside of both child pointers
    let mut inside_left = false;
    let mut inside_right = false;

    // Prepare an Intersections to collect filtered intersections
    let mut res_is = Intersections::new();

    for i in is.intersections.iter() {
        // If i.what is part of the left child, which_hit is true (for left)
        let which_hit = op.csg_left().includes(i.what);

        if intersection_allowed(op, which_hit, inside_left, inside_right) {
            res_is.intersections.push(i.clone());
        }

        // Depending on which object was hit, toggle either inside_left or
        // inside_right
        if which_hit {
            inside_left = !inside_left;
        } else {
            inside_right = !inside_right;
        }
    }
   
    // Return the filtered intersections
    res_is
}

#[test]
fn evaluating_rule_for_a_csg_operation() {
    // Dummy shapes for union CSG.
    let s1 = Shape::empty();
    let s2 = Shape::empty();

    let args: [(Shape, bool, bool, bool); 24] = [
        // Unions
        (Shape::csg_union(s1.clone(),s2.clone()),true,true,true),
        (Shape::csg_union(s1.clone(),s2.clone()),true,true,false),
        (Shape::csg_union(s1.clone(),s2.clone()),true,false,true),
        (Shape::csg_union(s1.clone(),s2.clone()),true,false,false),
        (Shape::csg_union(s1.clone(),s2.clone()),false,true,true),
        (Shape::csg_union(s1.clone(),s2.clone()),false,true,false),
        (Shape::csg_union(s1.clone(),s2.clone()),false,false,true),
        (Shape::csg_union(s1.clone(),s2.clone()),false,false,false),

        // Intersections
        (Shape::csg_intersection(s1.clone(),s2.clone()), true,true,true),
        (Shape::csg_intersection(s1.clone(),s2.clone()), true,true,false),
        (Shape::csg_intersection(s1.clone(),s2.clone()), true,false,true),
        (Shape::csg_intersection(s1.clone(),s2.clone()), true,false,false),
        (Shape::csg_intersection(s1.clone(),s2.clone()), false,true,true),
        (Shape::csg_intersection(s1.clone(),s2.clone()), false,true,false),
        (Shape::csg_intersection(s1.clone(),s2.clone()), false,false,true),
        (Shape::csg_intersection(s1.clone(),s2.clone()), false,false,false),

        // Differences
        (Shape::csg_difference(s1.clone(),s2.clone()), true,true,true),
        (Shape::csg_difference(s1.clone(),s2.clone()), true,true,false),
        (Shape::csg_difference(s1.clone(),s2.clone()), true,false,true),
        (Shape::csg_difference(s1.clone(),s2.clone()), true,false,false),
        (Shape::csg_difference(s1.clone(),s2.clone()), false,true,true),
        (Shape::csg_difference(s1.clone(),s2.clone()), false,true,false),
        (Shape::csg_difference(s1.clone(),s2.clone()), false,false,true),
        (Shape::csg_difference(s1.clone(),s2.clone()), false,false,false),
    ];

    let results: [bool; 24] = [
        // Unions
        false,
        true,
        false,
        true,
        false,
        false,
        true,
        true,

        // Intersections
        true,
        false,
        true,
        false,
        true,
        true,
        false,
        false,

        // Differences
        false,
        true,
        false,
        true,
        true,
        true,
        false,
        false
    ];

    for (arg, res) in args.iter().zip(results.iter()) {
        assert_eq!(intersection_allowed(&arg.0, arg.1, arg.2, arg.3), *res);
    }
}
