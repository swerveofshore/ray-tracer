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
#[derive(Clone, Debug)]
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

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Cube {
    pub transform: Matrix4D,
    pub material: Material,
}

impl Cube {
    pub fn new(transform: Matrix4D, material: Material) -> Cube {
        Cube { transform, material }
    }

    /// Creates a "unit cube."
    ///
    /// In object space, I consider a "unit cube" to be a cube with sides of
    /// length 1 in each direction, with the center of the cube at the origin.
    pub fn unit() -> Cube {
        Cube {
            transform: Default::default(),
            material: Default::default()
        }
    }

    /// Gets the minimum and maximum intersection offsets along an axis.
    ///
    /// No particular axis is specified; this function takes a component from
    /// the `origin` and `direction` fields of a `Ray4D` (e.g. `origin.x` and
    /// `direction.x`) and returns where the `Ray4D` intersects planes on a cube
    /// for a single axis.
    ///
    /// The smaller `t` is first in the tuple, the larger `t` is second.
    ///
    /// Note that this calculation assumes that the current `Cube` is a unit
    /// cube centered at the object-space origin.
    fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
        let tmin_numerator = -1.0 - origin;
        let tmax_numerator =  1.0 - origin;

        // Make sure that the direction is non-zero. If it is, assign INFINITY.
        let (tmin, tmax) = if direction.abs() >= FEQ_EPSILON {
            (tmin_numerator / direction, tmax_numerator / direction) 
        } else {
            (tmin_numerator * std::f64::INFINITY,
             tmax_numerator * std::f64::INFINITY)
        };

        // If tmin is actually greater than tmax, return tmax first.
        if tmin > tmax {
            (tmax, tmin)
        } else {
            (tmin, tmax)
        }
    }
}

impl Shape for Cube {
    fn local_intersect(&self, r: Ray4D) -> Intersections {
        let (xtmin, xtmax) = Cube::check_axis(r.origin.x, r.direction.x);
        let (ytmin, ytmax) = Cube::check_axis(r.origin.y, r.direction.y);
        let (ztmin, ztmax) = Cube::check_axis(r.origin.z, r.direction.z);

        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);

        // This can never happen if the ray hits the cube/sphere.
        if tmin > tmax {
            return Intersections::new()
        }

        Intersections {
            intersections: vec![
                Intersection { t: tmin, what: self },
                Intersection { t: tmax, what: self }
            ]
        }
    }

    /// Gets the normal for a face of the cube.
    ///
    /// Effectively, a normal is chosen by taking the component of `p` with the
    /// largest absolute value.
    ///
    /// If `p = (0.5, 0.7, -0.99)` for example, then the normal vector would be
    /// `(0.0, 0.0, -1)`, as that vector points in the direction of the face
    /// along the negative `z` axis.
    fn local_normal_at(&self, p: Tuple4D) -> Tuple4D {
        let xa = p.x.abs();
        let ya = p.y.abs();
        let za = p.z.abs();

        let max_component = xa.max(ya).max(za);
        if max_component == xa {
            Tuple4D::vector(p.x, 0.0, 0.0)
        } else if max_component == ya {
            Tuple4D::vector(0.0, p.y, 0.0)
        } else {
            Tuple4D::vector(0.0, 0.0, p.z)
        }
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

impl ShapeDebug for Cube { }

/// A cylinder.
///
/// Cylinders are centered along the Y axis. Specifically, the Y axis cuts
/// through the two "origins" of the cylinder's "circles" in object space.
///
/// By default, cylinders are "open," meaning that the circles on their ends
/// are not solid. In that case, they look somewhat like a paper towel roll.
///
/// To close the cylinder on its ends, set the `closed` property to `true`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Cylinder {
    pub transform: Matrix4D,
    pub material: Material,

    /// The bottom of the cylinder along the Y axis.
    pub minimum: f64,
    
    /// The top of the cylinder along the Y axis.
    pub maximum: f64,

    /// Whether the cylinder is closed.
    pub closed: bool,
}

impl Cylinder {
    pub fn unit() -> Cylinder {
        Cylinder {
            transform: Default::default(),
            material: Default::default(),

            minimum: -1.0 * std::f64::INFINITY,
            maximum: std::f64::INFINITY, 
            closed: false,
        }
    }

    pub fn intersect_caps<'a>(&'a self, r: Ray4D, is: &mut Intersections<'a>) {
        // If not closed, or the ray doesn't point anywhere near y, ignore caps.
        if !self.closed || r.direction.y.abs() < FEQ_EPSILON {
            return;
        }

        // Check for an intersection with the lower end cap.
        let tl = (self.minimum - r.origin.y) / r.direction.y;
        if Self::check_cap(r, tl) {
            is.intersections.push(Intersection { t: tl, what: self });
        }

        // Check for an intersection with the upper end cap.
        let tu = (self.maximum - r.origin.y) / r.direction.y;
        if Self::check_cap(r, tu) {
            is.intersections.push(Intersection { t: tu, what: self });
        }
    }

    /// Checks to see if the intersection `t` is within a radius of 1 (the
    /// assumed radius of a cylinder) from the Y axis.
    fn check_cap(r: Ray4D, t: f64) -> bool {
        let x = r.origin.x + t * r.direction.x;
        let z = r.origin.z + t * r.direction.z;

        (x.powi(2) + z.powi(2)) <= 1.0
    }
}

impl Shape for Cylinder {
    fn local_intersect(&self, r: Ray4D) -> Intersections {
        let a = r.direction.x.powi(2) + r.direction.z.powi(2); 
        
        // If r is parallel to the Y axis, check caps (if present) and leave.
        if a < FEQ_EPSILON {
            let mut is = Intersections::new();
            self.intersect_caps(r, &mut is);
            return is;
        }

        let b = 2.0f64 * r.origin.x * r.direction.x
              + 2.0f64 * r.origin.z * r.direction.z;

        let c = r.origin.x.powi(2) + r.origin.z.powi(2) - 1.0;

        let disc = b.powi(2) - 4.0 * a * c;

        // Ray r does not intersect the cylinder
        if disc < 0.0 {
            return Intersections::new();
        }

        // Ray r intersects the cylinder at one or two points
        let mut t0 = (-b - (disc.sqrt())) / (2.0 * a);
        let mut t1 = (-b + (disc.sqrt())) / (2.0 * a);
        
        // Make sure that t0 is the lowest intersection location.
        if t0 > t1 {
            std::mem::swap(&mut t0, &mut t1);
        }

        let mut is = Intersections::new();

        // If the t0 intersection is within the cylinder's bounds, add it.
        let y0 = r.origin.y + t0 * r.direction.y;
        if self.minimum < y0 && y0 < self.maximum {
            is.intersections.push(
                Intersection { t: t0, what: self }
            );
        }

        // If the t1 intersection is within the cylinder's bounds, add it.
        let y1 = r.origin.y + t1 * r.direction.y;
        if self.minimum < y1 && y1 < self.maximum {
            is.intersections.push(
                Intersection { t: t1, what: self }
            );
        }

        // Check whether any intersections occur at the cylinder caps.
        self.intersect_caps(r, &mut is);

        // Return all cylinder intersections.
        is
    }

    fn local_normal_at(&self, p: Tuple4D) -> Tuple4D {
        // Calculate the square of the distance from the y axis.
        let dist = p.x.powi(2) + p.z.powi(2);

        // If on the top cap, return a normal pointing up.
        if dist < 1.0 && p.y >= self.maximum - FEQ_EPSILON {
            Tuple4D::vector(0.0, 1.0, 0.0)
        }
        // If on the bottom cap, return a normal pointing down.
        else if dist < 1.0 && p.y <= self.minimum + FEQ_EPSILON {
            Tuple4D::vector(0.0, -1.0, 0.0)
        }
        // If on the round surface, return a normal pointing outwards.
        else {
            Tuple4D::vector(p.x, 0.0, p.z)
        }
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

impl ShapeDebug for Cylinder { }

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
    let p = Plane::new();

    let n1 = p.local_normal_at(Tuple4D::point(0.0, 0.0, 0.0));
    let n2 = p.local_normal_at(Tuple4D::point(10.0, 0.0, -10.0));
    let n3 = p.local_normal_at(Tuple4D::point(-5.0, 0.0, 150.0));

    assert_eq!(n1, Tuple4D::vector(0.0, 1.0, 0.0));
    assert_eq!(n2, Tuple4D::vector(0.0, 1.0, 0.0));
    assert_eq!(n3, Tuple4D::vector(0.0, 1.0, 0.0));
}

#[test]
fn ray_intersecting_plane_from_above() {
    let p = Plane::new();
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
    let p = Plane::new();
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

#[test]
fn under_point_is_below_the_surface() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let mut s = Sphere::glassy();
    s.transform = Matrix4D::translation(0.0, 0.0, 1.0);

    let i = Intersection { t: 5.0, what: &s };
    let is = Intersections { intersections: vec![i] };
    let comps = IntersectionComputation::new(&r, &i, Some(&is));

    assert!(comps.under_point.z > FEQ_EPSILON / 2.0);
    assert!(comps.point.z < comps.under_point.z);
}

#[test]
fn schlick_approximation_under_total_internal_reflection() {
    let s = Sphere::glassy();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 2.0f64.sqrt() / 2.0),
        Tuple4D::vector(0.0, 1.0, 0.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: -(2.0f64.sqrt()) / 2.0, what: &s },
            Intersection { t: 2.0f64.sqrt() / 2.0, what: &s }
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[1], Some(&is)
    );

    assert_eq!(1.0, comps.schlick());
}

#[test]
fn schlick_approximation_at_perpendicular_view() {
    use crate::feq;

    let s = Sphere::glassy();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.0),
        Tuple4D::vector(0.0, 1.0, 0.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: -1.0, what: &s },
            Intersection { t:  1.0, what: &s }
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[1], Some(&is)
    );

    assert!(feq(0.04, comps.schlick()));
}

#[test]
fn schlick_approximation_with_small_angle_and_n2_gt_n1() {
    use crate::feq;

    let s = Sphere::glassy();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.99, -2.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: 1.8589, what: &s },
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[0], Some(&is)
    );

    assert!(feq(0.48873, comps.schlick()));
}

#[test]
fn a_ray_intersects_a_cube() {
    let c = Cube::unit(); 
    let rays = vec![
        Ray4D::new(
            Tuple4D::point(5.0, 0.5, 0.0), Tuple4D::vector(-1.0, 0.0, 0.0)
        ),

        Ray4D::new(
            Tuple4D::point(-5.0, 0.5, 0.0), Tuple4D::vector(1.0, 0.0, 0.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.5, 5.0, 0.0), Tuple4D::vector(0.0, -1.0, 0.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.5, -5.0, 0.0), Tuple4D::vector(0.0, 1.0, 0.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.5, 0.0, 5.0), Tuple4D::vector(0.0, 0.0, -1.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.5, 0.0, -5.0), Tuple4D::vector(0.0, 0.0, 1.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 0.5, 0.0), Tuple4D::vector(0.0, 0.0, 1.0),
        )
    ];

    let expected_ts = vec![
        ( 4.0, 6.0),
        ( 4.0, 6.0),
        ( 4.0, 6.0),
        ( 4.0, 6.0),
        ( 4.0, 6.0),
        ( 4.0, 6.0),
        (-1.0, 1.0)
    ];

    // Show that local intersections are correctly calculated for any face.
    for (r, (t1, t2)) in rays.into_iter().zip(expected_ts.into_iter()) {
        let is = c.local_intersect(r);
        assert_eq!(is.intersections.len(), 2);
        assert_eq!(is.intersections[0].t, t1);
        assert_eq!(is.intersections[1].t, t2);
    }
}

#[test]
fn a_ray_misses_a_cube() {
    let c = Cube::unit(); 
    let rays = vec![
        Ray4D::new(
            Tuple4D::point(-2.0, 0.0, 0.0),
            Tuple4D::vector(0.2673, 0.5345, 0.8018)
        ),

        Ray4D::new(
            Tuple4D::point(0.0, -2.0, 0.0),
            Tuple4D::vector(0.8018, 0.2673, 0.5345)
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 0.0, -2.0),
            Tuple4D::vector(0.5345, 0.8018, 0.2673)
        ),

        Ray4D::new(
            Tuple4D::point(2.0, 0.0, 2.0),
            Tuple4D::vector(0.0, 0.0, -1.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 2.0, 2.0),
            Tuple4D::vector(0.0, -1.0, 0.0)
        ),

        Ray4D::new(
            Tuple4D::point(2.0, 2.0, 0.0),
            Tuple4D::vector(-1.0, 0.0, 0.0)
        ),
    ];

    // Show that local intersections are correctly calculated for any face.
    for r in rays {
        let is = c.local_intersect(r);
        assert_eq!(is.intersections.len(), 0);
    }
}

#[test]
fn normal_on_cube() {
    let c = Cube::unit();
    let points = vec![
        Tuple4D::point( 1.0,  0.5, -0.8),
        Tuple4D::point(-1.0, -0.2,  0.9),
        Tuple4D::point(-0.4,  1.0, -0.1),
        Tuple4D::point( 0.3, -1.0, -0.7),
        Tuple4D::point(-0.6,  0.3,  1.0),
        Tuple4D::point( 0.4,  0.4, -1.0),
        Tuple4D::point( 1.0,  1.0,  1.0),
        Tuple4D::point(-1.0, -1.0, -1.0),
    ];

    let normals = vec![
        Tuple4D::vector( 1.0,  0.0,  0.0),
        Tuple4D::vector(-1.0,  0.0,  0.0),
        Tuple4D::vector( 0.0,  1.0,  0.0),
        Tuple4D::vector( 0.0, -1.0,  0.0),
        Tuple4D::vector( 0.0,  0.0,  1.0),
        Tuple4D::vector( 0.0,  0.0, -1.0),
        Tuple4D::vector( 1.0,  0.0,  0.0),
        Tuple4D::vector(-1.0,  0.0,  0.0),
    ];

    for (p, n) in points.into_iter().zip(normals.into_iter()) {
        assert_eq!(n, c.local_normal_at(p));
    }
}

#[test]
fn a_ray_misses_a_cylinder() {
    let c = Cylinder::unit(); 
    let rays = vec![
        Ray4D::new(
            Tuple4D::point(1.0, 0.0, 0.0),
            Tuple4D::vector(0.0, 1.0, 0.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 0.0, 0.0),
            Tuple4D::vector(0.0, 1.0, 0.0)
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 0.0, -5.0),
            Tuple4D::vector(1.0, 1.0, 1.0)
        ),
    ];

    // Show that local intersections are correctly calculated for any face.
    for r in rays {
        let is = c.local_intersect(r);
        assert_eq!(is.intersections.len(), 0);
    }
}

#[test]
fn a_ray_strikes_a_cylinder() {
    use crate::feq;

    let c = Cylinder::unit(); 
    let rays = vec![
        Ray4D::new(
            Tuple4D::point(1.0, 0.0, -5.0),
            Tuple4D::vector(0.0, 0.0, 1.0).normalize()
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 0.0, -5.0),
            Tuple4D::vector(0.0, 0.0, 1.0).normalize()
        ),

        Ray4D::new(
            Tuple4D::point(0.5, 0.0, -5.0),
            Tuple4D::vector(0.1, 1.0, 1.0).normalize()
        ),
    ];

    let ts = vec![
        (5.0, 5.0),
        (4.0, 6.0),
        (6.80798, 7.08872)
    ];

    // Show that local intersections are correctly calculated for any face.
    for (r, (t0, t1)) in rays.into_iter().zip(ts.into_iter()) {
        let is = c.local_intersect(r);

        assert_eq!(is.intersections.len(), 2);
        assert!(feq(is.intersections[0].t, t0));
        assert!(feq(is.intersections[1].t, t1));
    }
}

#[test]
fn a_ray_strikes_a_constrained_cylinder() {
    use crate::feq;

    let mut c = Cylinder::unit(); 
    c.minimum = 1.0;
    c.maximum = 2.0;

    let rays = vec![
        Ray4D::new(
            Tuple4D::point(0.0, 1.5, 0.0),
            Tuple4D::vector(0.1, 1.0, 0.0).normalize()
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 3.0, -5.0),
            Tuple4D::vector(0.0, 0.0, 1.0).normalize()
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 0.0, -5.0),
            Tuple4D::vector(0.0, 0.0, 1.0).normalize()
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 2.0, -5.0),
            Tuple4D::vector(0.0, 0.0, 1.0).normalize()
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 1.0, -5.0),
            Tuple4D::vector(0.0, 0.0, 1.0).normalize()
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 1.5, -2.0),
            Tuple4D::vector(0.0, 0.0, 1.0).normalize()
        ),
    ];

    let counts = vec![ 0, 0, 0, 0, 0, 2 ];

    // Show that local intersections are correctly calculated for any face.
    for (r, co) in rays.into_iter().zip(counts.into_iter()) {
        let is = c.local_intersect(r);
        assert_eq!(is.intersections.len(), co);
    }
}

#[test]
fn a_ray_intersects_caps_of_a_closed_cylinder() {
    let mut c = Cylinder::unit();
    c.minimum = 1.0;
    c.maximum = 2.0;
    c.closed = true;

    let rays = vec![
        Ray4D::new(
            Tuple4D::point(0.0, 3.0, 0.0),
            Tuple4D::vector(0.0, -1.0, 0.0).normalize(),
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 3.0, -2.0),
            Tuple4D::vector(0.0, -1.0, 2.0).normalize(),
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 4.0, -2.0),
            Tuple4D::vector(0.0, -1.0, 1.0).normalize(),
        ),

        Ray4D::new(
            Tuple4D::point(0.0, 0.0, -2.0),
            Tuple4D::vector(0.0, 1.0, 2.0).normalize(),
        ),

        Ray4D::new(
            Tuple4D::point(0.0, -1.0, -2.0),
            Tuple4D::vector(0.0, 1.0, 1.0).normalize(),
        )
    ];

    let counts = vec![ 2, 2, 2, 2, 2 ];

    // Show that local intersections are correctly calculated for any face.
    for (r, co) in rays.into_iter().zip(counts.into_iter()) {
        let is = c.local_intersect(r);
        assert_eq!(is.intersections.len(), co);
    }
}
