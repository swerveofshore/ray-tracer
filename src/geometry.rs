use std::rc::Rc;
use std::cell::RefCell;

use crate::FEQ_EPSILON;
use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::ray::Ray4D;
use crate::light::Material;

/// An intersection.
///
/// This structure assumes that some ray produced an intersection. Parameter `t`
/// is analogous to `t` for a ray (the offset from the ray origin).
///
/// The `what` parameter is a reference to an `IntersectableDebug` trait object.
/// A trait object is used because there can be many unique objects (sphere,
/// cube, plane, etc.) which are intersectable.
#[derive(Clone, Debug)]
pub struct Intersection {
    pub t: f64,
    pub what: Rc<RefCell<dyn ShapeDebug>>,
}

/// Implements partial equality on an Intersection.
///
/// Two Intersection structures are equal if the offsets `t` of the
/// intersections are equivalent, and if the underlying *pointers* of the
/// intersections are equivalent.
impl PartialEq for Intersection {
    fn eq(&self, other: &Intersection) -> bool {
        self.t == other.t && Rc::ptr_eq(&self.what, &other.what)
    }
}

/// Implements total equality on an Intersection. Empty implementation.
impl Eq for Intersection { }

/// A collection of intersections.
///
/// Mostly a wrapper for a vector of `Intersection` objects. See the
/// `Intersection` documentation for more information.
pub struct Intersections {
    pub intersections: Vec<Intersection>,
}

impl Intersections {
    /// Creates a new list of intersections.
    pub fn new() -> Intersections {
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
    pub fn hit<'a>(&'a mut self) -> Option<&'a Intersection> {
        self.intersections.retain(|i| i.t.is_finite());
        self.sort();

        for i in self.intersections.iter() {
            if i.t >= 0.0 {
                return Some(i);
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

/// A combination of the Intersectable and the std::fmt::Debug trait.
pub trait ShapeDebug : Shape + std::fmt::Debug { }

/// Includes properties that all intersectable objects should have.
///
/// All intersectable objects should be able to be intersected by a ray.
/// Similarly, all intersectable objects should have normals on their surface.
pub trait Shape {
    fn intersect(&self, ray: Ray4D) -> Intersections;

    fn normal(&self, at: Tuple4D) -> Tuple4D;

    fn material(&self) -> &Material;
    fn material_mut(&mut self) -> &mut Material;

    fn transform(&self) -> &Matrix4D;
    fn transform_mut(&mut self) -> &mut Matrix4D;
}

/// A sphere.
///
/// A sphere is defined by its position, radius and transform. The `pos` field 
/// defines where a sphere is in *object* space, and the `transform` field moves
/// the sphere to *world* space.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Sphere {
    pub pos: Tuple4D,
    pub radius: f64,
    pub transform: Matrix4D,
    pub material: Material,
}

impl Sphere {
    /// Creates a new sphere with `radius` at `pos`.
    ///
    /// If `pos` is a vector, it is automatically converted to a point.
    pub fn new(mut pos: Tuple4D, radius: f64) -> Sphere {
        if !pos.is_point() {
            pos.w = 1.0;
        }

        Sphere {
            pos,
            radius,
            transform: Matrix4D::identity(),
            material: Default::default()
        }
    }

    /// Creates a unit sphere with default transform and material.
    pub fn unit() -> Sphere {
        Sphere {
            pos: Tuple4D::point(0.0, 0.0, 0.0),
            radius: 1.0,
            transform: Matrix4D::identity(),
            material: Default::default()
        }
    }
}

impl Shape for Sphere {
    /// Checks whether a ray intersects a sphere.
    ///
    /// Can return two values: either an empty vector (in the case of zero
    /// intersections), or a two-element vector (in the case of one or more
    /// intersections).
    ///
    /// If the sphere is only intersected at one point, the elements of the
    /// size-two vector should be equal.
    fn intersect(&self, ray: Ray4D) -> Intersections {
        let inverse_transform = self.transform.inverse().expect(
            "Transformation matrix on sphere should be invertible."
        );

        let transformed_ray = ray.transform(inverse_transform);
        let sphere_to_ray = transformed_ray.origin - self.pos; 

        let a = transformed_ray.direction.dot(&transformed_ray.direction);
        let b = 2.0 * transformed_ray.direction.dot(&sphere_to_ray);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b.powf(2.0) - (4.0 * a * c);

        if discriminant < 0.0 {
            return Intersections::new()
        }

        let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b + discriminant.sqrt()) / (2.0 * a);
        let i1 = Intersection { t: t1, what: Rc::new(RefCell::new(*self)) };
        let i2 = Intersection { t: t2, what: Rc::new(RefCell::new(*self)) };

        Intersections { intersections: vec![i1, i2] }
    }

    /// Returns the normal vector at point `at`.
    ///
    /// Note that a sphere specifies a transformation matrix. This matrix
    /// specifies the location of the sphere in world space, as well as any
    /// augmentations (scale, shearing, etc.) applied to the sphere.
    ///
    /// To compute the normal, `at` is first converted to object space.
    /// Afterwards, the normal is converted by taking a vector which points from
    /// the sphere's origin. Then, the vector is converted back to world space.
    ///
    /// The transposition of the transform inverse is taken to preserve the
    /// angle between the normal and the surface.
    fn normal(&self, at: Tuple4D) -> Tuple4D {
        let trans_inv = self.transform.inverse().expect(
            "Transformation matrix on sphere should be invertible."
        );

        let object_point = trans_inv * at;
        let object_normal = object_point - self.pos;
        let mut world_normal = trans_inv.transposition() * object_normal;
        world_normal.w = 0.0;

        world_normal.normalize()
    }

    /// Returns a reference to the material associated with a sphere.
    fn material(&self) -> &Material {
        &self.material
    }

    /// Returns a mutable reference to the material associated with a sphere.
    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }

    /// Returns a reference to the transform associated with a sphere.
    fn transform(&self) -> &Matrix4D {
        &self.transform
    }

    /// Returns a mutable reference to the transform associated with a sphere.
    fn transform_mut(&mut self) -> &mut Matrix4D {
        &mut self.transform
    }
}

/// Empty trait implementation for Sphere.
impl ShapeDebug for Sphere { }

pub struct Plane {
    pub normal: Tuple4D,
    pub material: Material,
    pub transform: Matrix4D,
}

impl Shape for Plane {
    fn normal(&self, _at: Tuple4D) -> Tuple4D {
        self.normal
    }

    fn intersect(&self, ray: Ray4D) -> Intersections {
        // TODO implement
        Intersections::new()
    }

    fn material(&self) -> &Material {
        &self.material
    }

    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }

    fn transform(&self) -> &Matrix4D {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Matrix4D {
        &mut self.transform
    }
}

/// A record for computations associated with an `Intersection`.
///
/// Mostly a superset of an `Intersection`.
#[derive(Clone, Debug)]
pub struct IntersectionComputation {
    /// The "time" of the ray intersection.
    pub t: f64,

    /// The object being intersected.
    pub obj: Rc<RefCell<dyn ShapeDebug>>,

    /// The point where the intersection occurs.
    pub point: Tuple4D,

    /// An error-corrected point of the intersection.
    pub over_point: Tuple4D,

    /// The eye vector for the intersection.
    pub eyev: Tuple4D,

    /// The normal vector of the object being intersected.
    pub normalv: Tuple4D,

    /// Whether the intersection occurs within the object or not.
    pub inside: bool,
}

impl IntersectionComputation {
    /// Creates a new intersection computation, given a ray and intersection.
    pub fn new(r: &Ray4D, i: &Intersection) -> IntersectionComputation {
        let t = i.t;
        let obj = Rc::clone(&i.what);
        let point = r.position(t);
        let eyev = -r.direction;
        let mut normalv = obj.borrow().normal(point);
        let over_point = point + normalv * FEQ_EPSILON;

        let inside = if normalv.dot(&eyev) < 0.0 {
            normalv = -normalv;
            true
        } else {
            false
        };

        IntersectionComputation {
            t, obj,
            point, over_point, eyev, normalv,
            inside
        }
    }
}

#[test]
fn ray_pierces_sphere() {
    let r = Ray4D::new(Tuple4D::point(0.0, 0.0, -5.0),
                       Tuple4D::vector(0.0, 0.0, 1.0));
    let s = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let xs = s.intersect(r);

    assert_eq!(xs.intersections.len(), 2);
    assert_eq!(xs.intersections[0].t, 4.0);
    assert_eq!(xs.intersections[1].t, 6.0);
}

#[test]
fn ray_is_tangent_to_sphere() {
    let r = Ray4D::new(Tuple4D::point(0.0, 1.0, -5.0),
                       Tuple4D::vector(0.0, 0.0, 1.0));
    let s = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let xs = s.intersect(r);

    assert_eq!(xs.intersections.len(), 2);
    assert_eq!(xs.intersections[0].t, 5.0);
    assert_eq!(xs.intersections[1].t, 5.0);
}

#[test]
fn ray_is_inside_sphere() {
    let r = Ray4D::new(Tuple4D::point(0.0, 0.0, 0.0),
                       Tuple4D::vector(0.0, 0.0, 1.0));
    let s = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let xs = s.intersect(r);

    assert_eq!(xs.intersections.len(), 2);
    assert_eq!(xs.intersections[0].t, -1.0);
    assert_eq!(xs.intersections[1].t, 1.0);
}

#[test]
fn sphere_is_behind_ray() {
    let r = Ray4D::new(Tuple4D::point(0.0, 0.0, 5.0),
                       Tuple4D::vector(0.0, 0.0, 1.0));
    let s = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let xs = s.intersect(r);

    assert_eq!(xs.intersections.len(), 2);
    assert_eq!(xs.intersections[0].t, -6.0);
    assert_eq!(xs.intersections[1].t, -4.0);
}

#[test]
fn hit_with_all_positive() {
    let s  = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let i1 = Intersection { t: 1.0, what: Rc::new(RefCell::new(s)) }; 
    let i2 = Intersection { t: 2.0, what: Rc::new(RefCell::new(s)) };
    let mut is = Intersections { intersections: vec![i1.clone(), i2.clone()] };

    assert_eq!(is.hit().unwrap(), &i1);
}

#[test]
fn hit_with_some_negative() {
    let s  = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let i1 = Intersection { t: -1.0, what: Rc::new(RefCell::new(s)) }; 
    let i2 = Intersection { t: 1.0, what: Rc::new(RefCell::new(s)) };
    let mut is = Intersections { intersections: vec![i1.clone(), i2.clone()] };

    assert_eq!(is.hit().unwrap(), &i2);
}

#[test]
fn hit_with_all_negative() {
    let s  = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let i1 = Intersection { t: -2.0, what: Rc::new(RefCell::new(s)) }; 
    let i2 = Intersection { t: -1.0, what: Rc::new(RefCell::new(s)) };
    let mut is = Intersections { intersections: vec![i1.clone(), i2.clone()] };

    assert_eq!(is.hit(), None);
}

#[test]
fn hit_multiple() {
    let s  = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let i1 = Intersection { t: 5.0,  what: Rc::new(RefCell::new(s)) }; 
    let i2 = Intersection { t: 7.0,  what: Rc::new(RefCell::new(s)) };
    let i3 = Intersection { t: -3.0, what: Rc::new(RefCell::new(s)) };
    let i4 = Intersection { t: 2.0,  what: Rc::new(RefCell::new(s)) };
    let mut is = Intersections {
        intersections: vec![i1.clone(), i2.clone(), i3.clone(), i4.clone()]
    };

    assert_eq!(is.hit().unwrap(), &i4);
}

#[test]
fn ray_hits_scaled_sphere() {
    let r = Ray4D::new(Tuple4D::point(0.0, 0.0, -5.0),
                       Tuple4D::vector(0.0, 0.0, 1.0));
    let mut s = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    s.transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let is = s.intersect(r);

    assert_eq!(is.intersections.len(), 2);
    assert_eq!(is.intersections[0].t, 3.0);
    assert_eq!(is.intersections[1].t, 7.0);
}

#[test]
fn ray_missess_translated_sphere() {
    let r = Ray4D::new(Tuple4D::point(0.0, 0.0, -5.0),
                       Tuple4D::vector(0.0, 0.0, 1.0));
    let mut s = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    s.transform = Matrix4D::translation(5.0, 0.0, 0.0);

    let is = s.intersect(r);
    assert_eq!(is.intersections.len(), 0);
}

#[test]
fn normal_on_sphere_x() {
    let s = Sphere::unit();
    let n = s.normal(Tuple4D::point(1.0, 0.0, 0.0));

    assert_eq!(n, Tuple4D::vector(1.0, 0.0, 0.0));
}

#[test]
fn normal_on_sphere_y() {
    let s = Sphere::unit();
    let n = s.normal(Tuple4D::point(0.0, 1.0, 0.0));

    assert_eq!(n, Tuple4D::vector(0.0, 1.0, 0.0));
}

#[test]
fn normal_on_sphere_z() {
    let s = Sphere::unit();
    let n = s.normal(Tuple4D::point(0.0, 0.0, 1.0));

    assert_eq!(n, Tuple4D::vector(0.0, 0.0, 1.0));
}

#[test]
fn normal_on_sphere_nonaxial() {
    let s = Sphere::unit();
    let n = s.normal(
        Tuple4D::point(
            3.0f64.sqrt() / 3.0,
            3.0f64.sqrt() / 3.0,
            3.0f64.sqrt() / 3.0
        )
    );

    assert_eq!(n,
        Tuple4D::vector(
            3.0f64.sqrt() / 3.0,
            3.0f64.sqrt() / 3.0,
            3.0f64.sqrt() / 3.0
        )
    );
}

#[test]
fn normal_on_sphere_translated() {
    let mut s = Sphere::unit();
    s.transform = Matrix4D::translation(0.0, 1.0, 0.0);
    let n = s.normal(Tuple4D::point(0.0, 1.70711, -0.70711));
    
    assert_eq!(n, Tuple4D::vector(0.0, 0.70711, -0.70711));
}

#[test]
fn normal_on_sphere_transformed() {
    let mut s = Sphere::unit();
    s.transform = Matrix4D::scaling(1.0, 0.5, 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 5.0);
    let n = s.normal(
        Tuple4D::point(0.0, 2.0f64.sqrt() / 2.0, -(2.0f64.sqrt() / 2.0))
    );

    assert_eq!(n, Tuple4D::vector(0.0, 0.97014, -0.24254));
}

#[test]
fn precompute_intersection_state() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let shape = Sphere::unit();
    let i = Intersection { t: 4.0, what: Rc::new(RefCell::new(shape)) };

    let comps = IntersectionComputation::new(&r, &i);

    assert!(Rc::ptr_eq(&comps.obj, &i.what));
    assert_eq!(comps.t, i.t);
    assert_eq!(comps.point, Tuple4D::point(0.0, 0.0, -1.0));
    assert_eq!(comps.eyev, Tuple4D::vector(0.0, 0.0, -1.0));
    assert_eq!(comps.normalv, Tuple4D::vector(0.0, 0.0, -1.0));
}

#[test]
fn precompute_outside_intersection() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let shape = Sphere::unit();
    let i = Intersection { t: 4.0, what: Rc::new(RefCell::new(shape)) };

    let comps = IntersectionComputation::new(&r, &i);
    assert!(!comps.inside);
}

#[test]
fn precompute_inside_intersection() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -0.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let shape = Sphere::unit();
    let i = Intersection { t: 1.0, what: Rc::new(RefCell::new(shape)) };

    let comps = IntersectionComputation::new(&r, &i);

    assert!(comps.inside);
    assert!(Rc::ptr_eq(&comps.obj, &i.what));
    assert_eq!(comps.t, i.t);
    assert_eq!(comps.point, Tuple4D::point(0.0, 0.0, 1.0));
    assert_eq!(comps.eyev, Tuple4D::vector(0.0, 0.0, -1.0));
    assert_eq!(comps.normalv, Tuple4D::vector(0.0, 0.0, -1.0));
}

#[test]
fn hit_should_offset_point() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let mut shape = Sphere::unit();
    shape.transform = Matrix4D::translation(0.0, 0.0, 1.0);

    let i = Intersection { t: 5.0, what: Rc::new(RefCell::new(shape)) };
    let comps = IntersectionComputation::new(&r, &i);

    assert!(comps.over_point.z < -FEQ_EPSILON / 2.0);
    assert!(comps.point.z > comps.over_point.z);
}
