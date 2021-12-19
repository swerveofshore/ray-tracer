use std::rc::Weak;

use crate::consts::FEQ_EPSILON;
use crate::tuple::Tuple4D;
use crate::ray::Ray4D;
use crate::light::Material;
use crate::matrix::Matrix4D;
use crate::intersect::{ Intersection, Intersections };

#[derive(Debug)]
pub enum ShapeType {
    /// An empty shape which does nothing. Mostly for testing.
    Empty,

    /// A unit sphere with its center at the object-space origin.
    Sphere,

    /// A 1-by-1-by-1 cube with its center at the object-space origin.
    Cube,

    /// A cone. Minimum Y, maximum Y and whether the cone is closed.
    Cone(f64, f64, bool),

    /// A cylinder. Minimum Y, maximum Y and whether the cylinder is closed.
    Cylinder(f64, f64, bool),

    /// A plane with a normal vector which stretches indefinitely along X and Y.
    Plane(Tuple4D),

    /// A group of shapes. Can include other groups of shapes.
    Group(Vec<Shape>),
}

#[derive(Debug)]
pub struct Shape {
    ty: ShapeType, 
    parent: Weak<Shape>,

    transform: Matrix4D,
    material: Material,
}

impl Shape {
    /// Creates a unit sphere with identity transform and default material.
    pub fn sphere() -> Shape {
        Shape {
            ty: ShapeType::Sphere,
            parent: Weak::new(),
            transform: Default::default(),
            material: Default::default(),
        }
    }

    /// Creates a plane with a normal pointing up along the Y axis.
    pub fn plane() -> Shape {
        Shape {
            ty: ShapeType::Plane(Tuple4D::vector(0.0, 1.0, 0.0)),
            parent: Weak::new(),
            transform: Default::default(),
            material: Default::default(),
        }
    }

    /// Creates a unit cube with identity transform and default material.
    pub fn cube() -> Shape {
        Shape {
            ty: ShapeType::Cube,
            parent: Weak::new(),
            transform: Default::default(),
            material: Default::default(),
        }
    }

    /// Creates an infinitely long cylinder with no end caps.
    pub fn cylinder() -> Shape {
        Shape {
            ty: ShapeType::Cylinder(
                    -1.0 * std::f64::INFINITY,
                    std::f64::INFINITY,
                    false
                ),
            parent: Weak::new(),
            transform: Default::default(),
            material: Default::default(),
        }
    }

    /// Creates a bounded cylinder without caps.
    pub fn bounded_cylinder(minimum: f64, maximum: f64) -> Shape {
         Shape {
            ty: ShapeType::Cylinder(minimum, maximum, false),
            parent: Weak::new(),
            transform: Default::default(),
            material: Default::default(),
        }       
    }

    /// Creates a bounded cylinder with caps.
    pub fn capped_cylinder(minimum: f64, maximum: f64) -> Shape {
         Shape {
            ty: ShapeType::Cylinder(minimum, maximum, true),
            parent: Weak::new(),
            transform: Default::default(),
            material: Default::default(),
        }       
    }

    /// Creates an infinitely long double-napped cone with no end caps.
    pub fn cone() -> Shape {
        Shape {
            ty: ShapeType::Cone(
                    -1.0 * std::f64::INFINITY,
                    std::f64::INFINITY,
                    false
                ),
            parent: Weak::new(),
            transform: Default::default(),
            material: Default::default(),
        }
    }

    pub fn transform(&self) -> &Matrix4D {
        &self.transform
    }

    pub fn transform_mut(&mut self) -> &mut Matrix4D {
        &mut self.transform
    }

    pub fn material(&self) -> &Material {
        &self.material
    }

    pub fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }

    pub fn local_intersect(&self, ray: &Ray4D) -> Intersections {
        match self.ty {
            ShapeType::Empty => Intersections::new(),
            ShapeType::Sphere => self.intersect_sphere(ray),
            ShapeType::Plane(_) => self.intersect_plane(ray),
            ShapeType::Cube => self.intersect_cube(ray),
            ShapeType::Cylinder(_, _, _) => self.intersect_cylinder(ray),
            ShapeType::Cone(_, _, _) => self.intersect_cone(ray),
            _ => Intersections::new()
        }
    }

    pub fn local_normal_at(&self, at: &Tuple4D) -> Tuple4D {
        match self.ty {
            ShapeType::Empty => Tuple4D { w: 0.0, ..*at },
            ShapeType::Sphere => self.normal_at_sphere(at),
            ShapeType::Plane(_) => self.normal_at_plane(at),
            ShapeType::Cube => self.normal_at_cube(at),
            ShapeType::Cylinder(_, _, _) => self.normal_at_cylinder(at),
            ShapeType::Cone(_, _, _) => self.normal_at_cone(at),
            _ => Default::default(),
        }
    }

    /// Checks whether a ray intersects a Sphere.
    ///
    /// Can return two values: either an empty vector (in the case of zero
    /// intersections), or a two-element vector (in the case of one or more
    /// intersections).
    ///
    /// If the sphere is only intersected at one point, the elements of the
    /// size-two vector should be equal.
    fn intersect_sphere(&self, ray: &Ray4D) -> Intersections {
        // Panic if intersect_sphere is being called on a non-sphere.
        match self.ty {
            ShapeType::Sphere => (),
            _ => unreachable!(),
        }

        // Sphere is centered at object-space origin.
        let sphere_to_ray = ray.origin;

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

    /// Returns the normal vector at a point on a Sphere.
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
    fn normal_at_sphere(&self, at: &Tuple4D) -> Tuple4D {
        // Panic if normal_at_sphere is being called on a non-sphere.
        match self.ty {
            ShapeType::Sphere => (),
            _ => unreachable!(),
        }

        Tuple4D { w: 0.0, ..*at }
    }

    fn intersect_plane(&self, ray: &Ray4D) -> Intersections {
        // Extract plane normal, panic if this isn't a plane.
        let _normal = match self.ty {
            ShapeType::Plane(n) => n,
            _ => unreachable!(),
        };

        // In local space, without a Y component, the ray won't intersect.
        if ray.direction.y.abs() <= FEQ_EPSILON {
            return Intersections::new();
        }

        // Calculate the offset of the ray with its Y component.
        let t = -ray.origin.y / ray.direction.y;
        let i = Intersection { t, what: self };

        Intersections { intersections: vec![i] }
    }

    fn normal_at_plane(&self, _at: &Tuple4D) -> Tuple4D {
        // Extract plane normal, panic if this isn't a plane.
        match self.ty {
            ShapeType::Plane(n) => n.clone(),
            _ => unreachable!(),
        }
    }

    fn intersect_cube(&self, ray: &Ray4D) -> Intersections {
        // Panic if this isn't a cube.
        match self.ty {
            ShapeType::Cube => (),
            _ => unreachable!(),
        }

        let (xtmin, xtmax)
            = Self::check_cube_axis(ray.origin.x, ray.direction.x);
        let (ytmin, ytmax)
            = Self::check_cube_axis(ray.origin.y, ray.direction.y);
        let (ztmin, ztmax)
            = Self::check_cube_axis(ray.origin.z, ray.direction.z);

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

    fn normal_at_cube(&self, p: &Tuple4D) -> Tuple4D {
        // Panic if this isn't a cube.
        match self.ty {
            ShapeType::Cube => (),
            _ => unreachable!(),
        }

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

    fn intersect_cylinder(&self, ray: &Ray4D) -> Intersections {
        let (minimum, maximum) = match self.ty {
            ShapeType::Cylinder(min, max, _) => (min, max),
            _ => unreachable!(),
        };

        let a = ray.direction.x.powi(2) + ray.direction.z.powi(2); 
        
        // If r is parallel to the Y axis, check caps (if present) and leave.
        if a < FEQ_EPSILON {
            let mut is = Intersections::new();
            self.intersect_cylinder_caps(ray, &mut is);
            return is;
        }

        let b = 2.0f64 * ray.origin.x * ray.direction.x
              + 2.0f64 * ray.origin.z * ray.direction.z;

        let c = ray.origin.x.powi(2) + ray.origin.z.powi(2) - 1.0;

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
        let y0 = ray.origin.y + t0 * ray.direction.y;
        if minimum < y0 && y0 < maximum {
            is.intersections.push(
                Intersection { t: t0, what: self }
            );
        }

        // If the t1 intersection is within the cylinder's bounds, add it.
        let y1 = ray.origin.y + t1 * ray.direction.y;
        if minimum < y1 && y1 < maximum {
            is.intersections.push(
                Intersection { t: t1, what: self }
            );
        }

        // Check whether any intersections occur at the cylinder caps.
        self.intersect_cylinder_caps(ray, &mut is);

        // Return all cylinder intersections.
        is
    }

    fn normal_at_cylinder(&self, at: &Tuple4D) -> Tuple4D {
        let (minimum, maximum) = match self.ty {
            ShapeType::Cylinder(min, max, _) => (min, max),
            _ => unreachable!(),
        };

        // Calculate the square of the distance from the y axis.
        let dist = at.x.powi(2) + at.z.powi(2);

        // If on the top cap, return a normal pointing up.
        if dist < 1.0 && at.y >= maximum - FEQ_EPSILON {
            Tuple4D::vector(0.0, 1.0, 0.0)
        }
        // If on the bottom cap, return a normal pointing down.
        else if dist < 1.0 && at.y <= minimum + FEQ_EPSILON {
            Tuple4D::vector(0.0, -1.0, 0.0)
        }
        // If on the round surface, return a normal pointing outwards.
        else {
            Tuple4D::vector(at.x, 0.0, at.z)
        }
    }

    fn intersect_cone(&self, ray: &Ray4D) -> Intersections {
        let (minimum, maximum) = match self.ty {
            ShapeType::Cone(min, max, _) => (min, max),
            _ => unreachable!(),
        };

        let a = ray.direction.x.powi(2)
              - ray.direction.y.powi(2)
              + ray.direction.z.powi(2); 
        
        let b = 2.0f64 * ray.origin.x * ray.direction.x
              - 2.0f64 * ray.origin.y * ray.direction.y
              + 2.0f64 * ray.origin.z * ray.direction.z;

        let c = ray.origin.x.powi(2)
              - ray.origin.y.powi(2)
              + ray.origin.z.powi(2);

        if a.abs() < FEQ_EPSILON {
            let mut is = Intersections::new();
            self.intersect_cone_caps(ray, &mut is);

            // Ray misses when both a and b are 0.
            if b.abs() < FEQ_EPSILON {
                return is;
            }
            // If only a is 0, then there's a single point of intersection.
            else {
                let t = -c / (2.0 * b);
                is.intersections.push(Intersection { t, what: self });
                return is;
            }
        }

        let disc = b.powi(2) - 4.0 * a * c;

        // Ray r does not intersect the cone 
        if disc < 0.0 {
            return Intersections::new();
        }

        // Ray r intersects the cone at one or two points
        let mut t0 = (-b - (disc.sqrt())) / (2.0 * a);
        let mut t1 = (-b + (disc.sqrt())) / (2.0 * a);
        
        // Make sure that t0 is the lowest intersection location.
        if t0 > t1 {
            std::mem::swap(&mut t0, &mut t1);
        }

        let mut is = Intersections::new();

        // If the t0 intersection is within the cone's bounds, add it.
        let y0 = ray.origin.y + t0 * ray.direction.y;
        if minimum < y0 && y0 < maximum {
            is.intersections.push(
                Intersection { t: t0, what: self }
            );
        }

        // If the t1 intersection is within the cone's bounds, add it.
        let y1 = ray.origin.y + t1 * ray.direction.y;
        if minimum < y1 && y1 < maximum {
            is.intersections.push(
                Intersection { t: t1, what: self }
            );
        }

        // Check whether any intersections occur at the cylinder caps.
        self.intersect_cone_caps(ray, &mut is);

        // Return all cylinder intersections.
        is
    }

    fn normal_at_cone(&self, at: &Tuple4D) -> Tuple4D {
        let (minimum, maximum) = match self.ty {
            ShapeType::Cone(min, max, _) => (min, max),
            _ => unreachable!(),
        };

        // Calculate the square of the distance from the y axis.
        let dist = at.x.powi(2) + at.z.powi(2);

        // If on the top cap, return a normal pointing up.
        if dist < 1.0 && at.y >= maximum - FEQ_EPSILON {
            Tuple4D::vector(0.0, 1.0, 0.0)
        }
        // If on the bottom cap, return a normal pointing down.
        else if dist < 1.0 && at.y <= minimum + FEQ_EPSILON {
            Tuple4D::vector(0.0, -1.0, 0.0)
        }
        // If on the round surface, return a normal pointing outwards.
        else {
            let mut y = dist.sqrt();
            if at.y > 0.0 {
                y = -y;
            }

            Tuple4D::vector(at.x, y, at.z)
        }
    }

    /// Gets the min and max intersection offsets along an axis of a cube.
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
    fn check_cube_axis(origin: f64, direction: f64) -> (f64, f64) {
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

    fn intersect_cylinder_caps<'a>(&'a self, ray: &Ray4D,
        is: &mut Intersections<'a>) {
        let (minimum, maximum, closed) = match self.ty {
            ShapeType::Cylinder(min, max, c) => (min, max, c),
            _ => unreachable!(),
        };

        // If not closed, or the ray doesn't point anywhere near y, ignore caps.
        if !closed || ray.direction.y.abs() < FEQ_EPSILON {
            return;
        }

        // Check for an intersection with the lower end cap.
        let tl = (minimum - ray.origin.y) / ray.direction.y;
        if Self::check_cylinder_cap(ray, tl) {
            is.intersections.push(Intersection { t: tl, what: self });
        }

        // Check for an intersection with the upper end cap.
        let tu = (maximum - ray.origin.y) / ray.direction.y;
        if Self::check_cylinder_cap(ray, tu) {
            is.intersections.push(Intersection { t: tu, what: self });
        }
    }

    /// Checks to see if the intersection `t` is within a radius of 1 (the
    /// assumed radius of a cylinder) from the Y axis.
    fn check_cylinder_cap(ray: &Ray4D, t: f64) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;

        (x.powi(2) + z.powi(2)) <= 1.0
    }

    fn intersect_cone_caps<'a>(&'a self, ray: &Ray4D,
        is: &mut Intersections<'a>) {
        let (minimum, maximum, closed) = match self.ty {
            ShapeType::Cone(min, max, c) => (min, max, c),
            _ => unreachable!(),
        };

        // If not closed, or the ray doesn't point anywhere near y, ignore caps.
        if closed || ray.direction.y.abs() < FEQ_EPSILON {
            return;
        }

        // Check for an intersection with the lower end cap.
        let tl = (minimum - ray.origin.y) / ray.direction.y;
        if Self::check_cone_cap(ray, tl, minimum) {
            is.intersections.push(Intersection { t: tl, what: self });
        }

        // Check for an intersection with the upper end cap.
        let tu = (maximum - ray.origin.y) / ray.direction.y;
        if Self::check_cone_cap(ray, tu, maximum) {
            is.intersections.push(Intersection { t: tu, what: self });
        }
    }

    fn check_cone_cap(ray: &Ray4D, t: f64, y: f64) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;

        x.powi(2) + z.powi(2) <= y.powi(2)
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
pub fn intersect(s: &Shape, r: Ray4D) -> Intersections {
    let inverse_transform = s.transform().inverse().expect(
        "Transformation matrix on shape should be invertible."
    );

    let transformed_ray = r.transform(inverse_transform);
    s.local_intersect(&transformed_ray)
}

pub fn normal_at(s: &Shape, p: Tuple4D) -> Tuple4D {
    let trans_inv = s.transform.inverse().expect(
        "Shape transform should be invertible."
    );

    let local_point = trans_inv * p;
    let local_normal = s.local_normal_at(&local_point);
    let mut world_normal = trans_inv.transposition() * local_normal;
    world_normal.w = 0.0;

    world_normal.normalize()
}