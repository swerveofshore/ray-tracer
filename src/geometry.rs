use crate::consts::FEQ_EPSILON;
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

/// Implements total equality on an Intersection. Empty implementation.
impl<'a> Eq for Intersection<'a> { }

/// A collection of intersections.
///
/// Mostly a wrapper for a vector of `Intersection` objects. See the
/// `Intersection` documentation for more information.
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

    pub fn hit_imm(&self) -> Option<Intersection<'a>> {
        let mut real_intersections: Vec<&Intersection<'a>> =
            self.intersections.iter().filter(|i| i.t.is_finite()).collect();
        real_intersections.sort_by(|a, b|
            a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal)
        );

        for i in real_intersections {
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

/// Includes properties that all intersectable objects should have.
///
/// All intersectable objects should be able to be intersected by a ray.
/// Similarly, all intersectable objects should have normals on their surface.
pub trait Shape {
    fn local_intersect(&self, ray: Ray4D) -> Intersections;
    fn local_normal_at(&self, at: Tuple4D) -> Tuple4D;

    fn material(&self) -> &Material;
    fn material_mut(&mut self) -> &mut Material;

    fn transform(&self) -> &Matrix4D;
    fn transform_mut(&mut self) -> &mut Matrix4D;
}

/// A combination of the Intersectable and the std::fmt::Debug trait.
pub trait ShapeDebug : Shape + std::fmt::Debug { }

/// A sphere.
///
/// A sphere is defined by its position, radius and transform. The `pos` field 
/// defines where a sphere is in *object* space, and the `transform` field moves
/// the sphere to *world* space.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Sphere {
    pub radius: f64,
    pub pos: Tuple4D,

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
            radius,
            pos,
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    /// Creates a unit sphere with default transform and material.
    pub fn unit() -> Sphere {
        Sphere {
            pos: Tuple4D::point(0.0, 0.0, 0.0),
            radius: 1.0,
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    /// Creates a unit sphere with a glassy material.
    pub fn glassy() -> Sphere {
        Sphere {
            material: Material {
                transparency: 1.0,
                refractive_index: 1.5,
                ..Default::default()
            },

            ..Sphere::unit()
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
    fn local_intersect(&self, ray: Ray4D) -> Intersections {
        let sphere_to_ray = ray.origin - self.pos; 

        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * ray.direction.dot(&sphere_to_ray);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b.powf(2.0) - (4.0 * a * c);

        if discriminant < 0.0 {
            return Intersections::new()
        }

        let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b + discriminant.sqrt()) / (2.0 * a);
        let i1 = Intersection { t: t1, what: self };
        let i2 = Intersection { t: t2, what: self };

        Intersections { intersections: vec![i1, i2] }
    }

    /// Returns the normal vector at point `at`.
    ///
    /// Note that the following calculations occur in object space. This
    /// function assumes that `at` is local to the sphere. Function `normal_at`
    /// does a conversion from world space to local space before calling this
    /// function, for example.
    ///
    /// A normal vector on the surface of a sphere can be found by subtracting
    /// the sphere's position (origin) from a point on the surface.
    ///
    /// Point `at` should be a point; the value returned from this function is
    /// a vector.
    fn local_normal_at(&self, at: Tuple4D) -> Tuple4D {
        let mut sphere_normal = at - self.pos;
        sphere_normal.w = 0.0;

        sphere_normal
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

/// A plane.
///
/// Defines a plane which stretches indefinitely. The orientation of the plane
/// is defined by a normal vector.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Plane {
    pub normal: Tuple4D,

    pub material: Material,
    pub transform: Matrix4D,
}

impl Plane {
    pub fn new() -> Plane {
        Plane {
            normal: Tuple4D::vector(0.0, 1.0, 0.0),
            material: Default::default(),
            transform: Default::default(),
        }
    }
}

impl Shape for Plane {
    fn local_intersect(&self, ray: Ray4D) -> Intersections {
        // In local space, plane stretches across XZ axis; if the ray doesn't
        // point in the Y direction whatsoever, it won't intersect the plane.
        if ray.direction.y.abs() <= FEQ_EPSILON {
            return Intersections::new();
        }
        
        // In local space, we can deduce the offset of the ray by checking
        // the Y component. If a ray points down, it intersects with positive t.
        let t = -ray.origin.y / ray.direction.y;
        let i = Intersection { t, what: self };

        Intersections { intersections: vec![i] }
    }

    fn local_normal_at(&self, _at: Tuple4D) -> Tuple4D {
        self.normal
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

impl ShapeDebug for Plane { }

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

    /// An error-corrected point of the intersection.
    pub over_point: Tuple4D,

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
        let over_point = point + normalv * FEQ_EPSILON;

        let inside = if normalv.dot(&eyev) < 0.0 {
            normalv = -normalv;
            true
        } else {
            false
        };

        let reflectv = r.direction.reflect(&normalv);
        let (n1, n2) = if let Some(xs) = is {
            Self::refraction_indices(hit, xs)
        } else {
            (1.0, 1.0)
        };

        IntersectionComputation {
            t, obj,
            point, over_point,
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
}

/// An empty shape for testing.
///
/// Consider this shape a "null" shape which doesn't define any points or
/// surfaces. Literally it is nothing in space, yet it still defines a transform
/// and material (and can save rays).
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct EmptyShape {
    pub transform: Matrix4D,
    pub material: Material,
}

impl EmptyShape {
    #[allow(unused)] // Allowed because EmptyShape is mostly for testing.
    pub(crate) fn new() -> EmptyShape {
        EmptyShape {
            transform: Default::default(),
            material: Default::default(),
        }
    }
}

impl Shape for EmptyShape {
    fn local_intersect(&self, _r: Ray4D) -> Intersections {
        // Empty shape can never be intersected.
        Intersections::new()
    }

    fn local_normal_at(&self, at: Tuple4D) -> Tuple4D {
        // Assume that an empty object is a sphere at the origin.
        let mut normal = at;
        normal.w = 0.0;

        normal
    }

    fn transform(&self) -> &Matrix4D {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Matrix4D {
        &mut self.transform
    }

    fn material(&self) -> &Material {
        &self.material
    }

    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }
}

impl ShapeDebug for EmptyShape { }

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

/// Calculate the normal vector for a `Shape`.
///
/// Each `Shape` implements `local_normal_at`, which calculates a shape's normal
/// in *object* space. This function uses the `transform` associated with each
/// `Shape`, and converts the object-space normal to a world-space normal.
///
/// Subsequently, each world-space normal is properly orthogonal to the `Shape`
/// in world space. Moreover, the normal is normalized before being returned.
pub fn normal_at(s: & dyn ShapeDebug, p: Tuple4D) -> Tuple4D {
    let trans_inv = s.transform().inverse().expect(
        "Shape transform should be invertible."
    );

    let local_point = trans_inv * p;
    let local_normal = s.local_normal_at(local_point);
    let mut world_normal = trans_inv.transposition() * local_normal;
    world_normal.w = 0.0;

    world_normal.normalize()
}

#[test]
fn compute_normal_on_translated_sphere() {
    let mut s = Sphere::unit();
    s.transform = Matrix4D::translation(0.0, 1.0, 0.0);

    let p = Tuple4D::point(0.0, 1.70711, -0.70711);
    let n = normal_at(&s, p);

    assert_eq!(n, Tuple4D::vector(0.0, 0.70711, -0.70711));
}

#[test]
fn compute_normal_on_transformed_sphere() {
    let mut s = Sphere::unit();
    s.transform = Matrix4D::scaling(1.0, 0.5, 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 5.0);
    
    let p = Tuple4D::point(0.0, 2.0f64.sqrt() / 2.0, -(2.0f64.sqrt()) / 2.0);
    let n = normal_at(&s, p);

    assert_eq!(n, Tuple4D::vector(0.0, 0.97014, -0.24254));
}

#[test]
fn normal_on_plane() {
    let mut p = Plane::new();

    let n1 = p.local_normal_at(Tuple4D::point(0.0, 0.0, 0.0));
    let n2 = p.local_normal_at(Tuple4D::point(10.0, 0.0, -10.0));
    let n3 = p.local_normal_at(Tuple4D::point(-5.0, 0.0, 150.0));

    assert_eq!(n1, Tuple4D::vector(0.0, 1.0, 0.0));
    assert_eq!(n2, Tuple4D::vector(0.0, 1.0, 0.0));
    assert_eq!(n3, Tuple4D::vector(0.0, 1.0, 0.0));
}

#[test]
fn ray_intersecting_plane_from_above() {
    let mut p = Plane::new();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 1.0, 0.0),
        Tuple4D::vector(0.0, -1.0, 0.0)
    );

    let is = p.local_intersect(r);

    assert_eq!(is.intersections.len(), 1);
    assert_eq!(is.intersections[0].t, 1.0);
}

#[test]
fn ray_intersecting_plane_from_below() {
    let mut p = Plane::new();
    let r = Ray4D::new(
        Tuple4D::point(0.0, -1.0, 0.0),
        Tuple4D::vector(0.0, 1.0, 0.0)
    );

    let is = p.local_intersect(r);

    assert_eq!(is.intersections.len(), 1);
    assert_eq!(is.intersections[0].t, 1.0);
}

/*
#[test]
fn ray_is_tangent_to_sphere() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 1.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );
    let s = Sphere::new(Tuple4D::point(0.0, 0.0, 0.0), 1.0);
    let xs = intersect(&*(s.borrow()), r);

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
*/

#[test]
fn find_refraction_indices_from_intersections() {
    let mut a = Sphere::glassy();
    a.transform = Matrix4D::scaling(2.0, 2.0, 2.0);
    a.material.refractive_index = 1.5;

    let mut b = Sphere::glassy();
    b.transform = Matrix4D::translation(0.0, 0.0, -0.25);
    b.material.refractive_index = 2.0;

    let mut c = Sphere::glassy();
    c.transform = Matrix4D::translation(0.0, 0.0, 0.25);
    c.material.refractive_index = 2.5;

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -4.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: 2.0,  what: &a },
            Intersection { t: 2.75, what: &b },
            Intersection { t: 3.25, what: &c },
            Intersection { t: 4.75, what: &b },
            Intersection { t: 5.25, what: &c },
            Intersection { t: 6.0,  what: &a },
        ],
    };

    let expected_refractions = vec![
        (1.0, 1.5),
        (1.5, 2.0),
        (2.0, 2.5),
        (2.5, 2.5),
        (2.5, 1.5),
        (1.5, 1.0)
    ];

    for (j, i) in is.intersections.iter().enumerate() {
        let comps = IntersectionComputation::new(&r, &i, Some(&is));
        assert_eq!(expected_refractions[j], (comps.n1, comps.n2)); 
    }
}
