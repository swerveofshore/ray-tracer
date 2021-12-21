use std::rc::{ Rc, Weak };
use std::cell::{ Ref, RefCell };

use crate::consts::FEQ_EPSILON;
use crate::tuple::Tuple4D;
use crate::ray::Ray4D;
use crate::light::Material;
use crate::matrix::Matrix4D;
use crate::intersect::{ Intersection, Intersections };

#[derive(Debug)]
pub enum ShapeNodeType {
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
    Group(Vec<Rc<RefCell<ShapeNode>>>),
}

#[derive(Debug)]
pub struct ShapeNode {
    pub ty: ShapeNodeType, 
    parent: Weak<RefCell<ShapeNode>>,

    pub transform: Matrix4D,
    pub material: Material,
}

impl ShapeNode {
    /// Creates an empty shape which does nothing. Mostly for testing.
    pub fn empty() -> ShapeNode {
        ShapeNode {
            ty: ShapeNodeType::Empty,
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    /// Creates a unit sphere with identity transform and default material.
    pub fn sphere() -> ShapeNode {
        ShapeNode {
            ty: ShapeNodeType::Sphere,
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    /// Creates a plane with a normal pointing up along the Y axis.
    pub fn plane() -> ShapeNode {
        ShapeNode {
            ty: ShapeNodeType::Plane(Tuple4D::vector(0.0, 1.0, 0.0)),
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    /// Creates a unit cube with identity transform and default material.
    pub fn cube() -> ShapeNode {
        ShapeNode {
            ty: ShapeNodeType::Cube,
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    /// Creates an infinitely long cylinder with no end caps.
    pub fn cylinder() -> ShapeNode {
        ShapeNode {
            ty: ShapeNodeType::Cylinder(
                    -1.0 * std::f64::INFINITY,
                    std::f64::INFINITY,
                    false
                ),
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    /// Creates a bounded cylinder without caps.
    pub fn bounded_cylinder(minimum: f64, maximum: f64) -> ShapeNode {
         ShapeNode {
            ty: ShapeNodeType::Cylinder(minimum, maximum, false),
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }       
    }

    /// Creates a bounded cylinder with caps.
    pub fn capped_cylinder(minimum: f64, maximum: f64) -> ShapeNode {
         ShapeNode {
            ty: ShapeNodeType::Cylinder(minimum, maximum, true),
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }       
    }

    /// Creates an infinitely long double-napped cone with no end caps.
    pub fn cone() -> ShapeNode {
        ShapeNode {
            ty: ShapeNodeType::Cone(
                    -1.0 * std::f64::INFINITY,
                    std::f64::INFINITY,
                    false
                ),
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
        }
    }

    pub fn group() -> ShapeNode {
        ShapeNode {
            ty: ShapeNodeType::Group(Vec::new()),
            parent: Weak::new(),
            transform: Matrix4D::identity(),
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

    pub fn world_to_object(&self, mut point: Tuple4D) -> Tuple4D {
        // If the current object has a parent, convert the parent point to
        // object space first.
        if let Some(parent) = self.parent.upgrade() {
            point = parent.borrow().world_to_object(point);
        }

        let trans_inv = self.transform.inverse().expect(
            "ShapeNode transform should be invertible."
        );

        // Then, convert the leaf (youngest child) point to object space.
        trans_inv * point
    }

    pub fn normal_to_world(&self, mut normal: Tuple4D) -> Tuple4D {
        let trans_inv = self.transform.inverse().expect(
            "ShapeNode transform should be invertible."
        );

        normal = trans_inv.transposition() * normal;
        normal.w = 0.0;
        normal = normal.normalize();

        // If the current object has a parent, convert the normal to the
        // parent's world space.
        if let Some(parent) = self.parent.upgrade() {
            normal = parent.borrow().normal_to_world(normal);
        }

        normal
    }

    pub fn local_intersect(&self, ray: &Ray4D) -> Intersections {
        match self.ty {
            ShapeNodeType::Empty => Intersections::new(),
            ShapeNodeType::Sphere => self.intersect_sphere(ray),
            ShapeNodeType::Plane(_) => self.intersect_plane(ray),
            ShapeNodeType::Cube => self.intersect_cube(ray),
            ShapeNodeType::Cylinder(_, _, _) => self.intersect_cylinder(ray),
            ShapeNodeType::Cone(_, _, _) => self.intersect_cone(ray),
            ShapeNodeType::Group(_) => self.intersect_group(ray),
        }
    }

    pub fn local_normal_at(&self, at: &Tuple4D) -> Tuple4D {
        match self.ty {
            ShapeNodeType::Empty => Tuple4D { w: 0.0, ..*at },
            ShapeNodeType::Sphere => self.normal_at_sphere(at),
            ShapeNodeType::Plane(_) => self.normal_at_plane(at),
            ShapeNodeType::Cube => self.normal_at_cube(at),
            ShapeNodeType::Cylinder(_, _, _) => self.normal_at_cylinder(at),
            ShapeNodeType::Cone(_, _, _) => self.normal_at_cone(at),
            ShapeNodeType::Group(_) => panic!(
                "Local normal calculations should never occur on groups."
            ),
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
            ShapeNodeType::Sphere => (),
            _ => unreachable!(),
        }

        // Sphere is centered at object-space origin.
        // Note that subtracting a point removes the 'w' part of the ray origin.
        let sphere_to_ray = ray.origin - Tuple4D::point(0.0, 0.0, 0.0);

        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * ray.direction.dot(&sphere_to_ray);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b.powi(2) - (4.0 * a * c);

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
            ShapeNodeType::Sphere => (),
            _ => unreachable!(),
        }

        Tuple4D { w: 0.0, ..*at }
    }

    fn intersect_plane(&self, ray: &Ray4D) -> Intersections {
        // Extract plane normal, panic if this isn't a plane.
        let _normal = match self.ty {
            ShapeNodeType::Plane(n) => n,
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
            ShapeNodeType::Plane(n) => n.clone(),
            _ => unreachable!(),
        }
    }

    fn intersect_cube(&self, ray: &Ray4D) -> Intersections {
        // Panic if this isn't a cube.
        match self.ty {
            ShapeNodeType::Cube => (),
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
            ShapeNodeType::Cube => (),
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
            ShapeNodeType::Cylinder(min, max, _) => (min, max),
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
            ShapeNodeType::Cylinder(min, max, _) => (min, max),
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
            ShapeNodeType::Cone(min, max, _) => (min, max),
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
            ShapeNodeType::Cone(min, max, _) => (min, max),
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

    fn intersect_group(&self, ray: &Ray4D) -> Intersections {
        let children = match self.ty {
            ShapeNodeType::Group(ref c) => c,
            _ => unreachable!(),
        };

        // If there are no children, return an empty list of intersections.
        if children.is_empty() {
            return Intersections::new()
        }

        let children_refs: Vec<_> = children.iter().map(
            // Leak the child reference--we are *certain* that it will persist.
            // As an aside, I'm not super familiar with what Ref::leak()
            // really does, so this may need to be substituted later.
            |child| Ref::leak(child.borrow())
        ).collect();

        // Otherwise, for each child, collect its intersections and aggregate.
        let mut all_intersections = Vec::new();
        for cr in children_refs.iter() {
            // Use regular intersect, otherwise transforms won't be applied
            let is = intersect(cr, *ray);
            all_intersections.push(is)
        }

        Intersections::aggregate(all_intersections)
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
            ShapeNodeType::Cylinder(min, max, c) => (min, max, c),
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
            ShapeNodeType::Cone(min, max, c) => (min, max, c),
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

pub fn add_child_to_group(group: Rc<RefCell<ShapeNode>>,
    child: Rc<RefCell<ShapeNode>>) {

    let mut gb = group.borrow_mut();
    let children = match gb.ty {
        ShapeNodeType::Group(ref mut c) => c,
        _ => panic!("Cannot add child to non-group shape."),
    };

    child.borrow_mut().parent = Rc::downgrade(&group);
    children.push(child);
}

/// Intersects a ray with a `ShapeNode`.
///
/// Each shape implements `local_intersect`, which calculates intersections of a
/// shape in *normal* space. This function uses the `transform` associated with
/// each `ShapeNode`, and converts the ray to object space before calculating
/// intersections.
///
/// Intersections are technically calculated in local space; these local-space
/// intersections are returned in an `Intersections` record.
pub fn intersect(s: &ShapeNode, r: Ray4D) -> Intersections {
    let inverse_transform = s.transform().inverse().expect(
        "Transformation matrix on shape should be invertible."
    );

    let transformed_ray = r.transform(inverse_transform);
    s.local_intersect(&transformed_ray)
}

pub fn normal_at(s: &ShapeNode, world_point: Tuple4D) -> Tuple4D {
    let local_point = s.world_to_object(world_point);
    let local_normal = s.local_normal_at(&local_point);
    s.normal_to_world(local_normal)
}

#[test]
fn compute_normal_on_translated_sphere() {
    let mut s = ShapeNode::sphere();
    s.transform = Matrix4D::translation(0.0, 1.0, 0.0);

    let p = Tuple4D::point(0.0, 1.70711, -0.70711);
    let n = normal_at(&s, p);

    assert_eq!(n, Tuple4D::vector(0.0, 0.70711, -0.70711));
}

#[test]
fn compute_normal_on_transformed_sphere() {
    let mut s = ShapeNode::sphere();
    s.transform = Matrix4D::scaling(1.0, 0.5, 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 5.0);
    
    let p = Tuple4D::point(0.0, 2.0f64.sqrt() / 2.0, -(2.0f64.sqrt()) / 2.0);
    let n = normal_at(&s, p);

    assert_eq!(n, Tuple4D::vector(0.0, 0.97014, -0.24254));
}

#[test]
fn normal_on_plane() {
    let p = ShapeNode::plane();

    let n1 = p.local_normal_at(&Tuple4D::point(0.0, 0.0, 0.0));
    let n2 = p.local_normal_at(&Tuple4D::point(10.0, 0.0, -10.0));
    let n3 = p.local_normal_at(&Tuple4D::point(-5.0, 0.0, 150.0));

    assert_eq!(n1, Tuple4D::vector(0.0, 1.0, 0.0));
    assert_eq!(n2, Tuple4D::vector(0.0, 1.0, 0.0));
    assert_eq!(n3, Tuple4D::vector(0.0, 1.0, 0.0));
}

#[test]
fn ray_intersecting_plane_from_above() {
    let p = ShapeNode::plane();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 1.0, 0.0),
        Tuple4D::vector(0.0, -1.0, 0.0)
    );

    let is = p.local_intersect(&r);

    assert_eq!(is.intersections.len(), 1);
    assert_eq!(is.intersections[0].t, 1.0);
}

#[test]
fn ray_intersecting_plane_from_below() {
    let p = ShapeNode::plane();
    let r = Ray4D::new(
        Tuple4D::point(0.0, -1.0, 0.0),
        Tuple4D::vector(0.0, 1.0, 0.0)
    );

    let is = p.local_intersect(&r);

    assert_eq!(is.intersections.len(), 1);
    assert_eq!(is.intersections[0].t, 1.0);
}

#[test]
fn ray_is_tangent_to_sphere() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 1.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );
    let s = ShapeNode::sphere();
    let xs = intersect(&s, r);

    assert_eq!(xs.intersections.len(), 2);
    assert_eq!(xs.intersections[0].t, 5.0);
    assert_eq!(xs.intersections[1].t, 5.0);
}

#[test]
fn ray_is_inside_sphere() {
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let s = ShapeNode::sphere();

    let xs = s.intersect_sphere(&r);
    assert_eq!(xs.intersections.len(), 2);
    assert_eq!(xs.intersections[0].t, -1.0);
    assert_eq!(xs.intersections[1].t, 1.0);
}

#[test]
fn sphere_is_behind_ray() {
    let r = Ray4D::new(Tuple4D::point(0.0, 0.0, 5.0),
                       Tuple4D::vector(0.0, 0.0, 1.0));
    let s = ShapeNode::sphere();
    let xs = intersect(&s, r);

    assert_eq!(xs.intersections.len(), 2);
    assert_eq!(xs.intersections[0].t, -6.0);
    assert_eq!(xs.intersections[1].t, -4.0);
}

#[test]
fn hit_with_all_positive() {
    let s  = ShapeNode::sphere();
    let i1 = Intersection { t: 1.0, what: &s }; 
    let i2 = Intersection { t: 2.0, what: &s };
    let mut is = Intersections { intersections: vec![i1, i2] };

    assert_eq!(is.hit().unwrap(), i1);
}

#[test]
fn hit_with_some_negative() {
    let s  = ShapeNode::sphere();
    let i1 = Intersection { t: -1.0, what: &s }; 
    let i2 = Intersection { t: 1.0, what: &s };
    let mut is = Intersections { intersections: vec![i1, i2] };

    assert_eq!(is.hit().unwrap(), i2);
}

#[test]
fn hit_with_all_negative() {
    let s  = ShapeNode::sphere();
    let i1 = Intersection { t: -2.0, what: &s };
    let i2 = Intersection { t: -1.0, what: &s };
    let mut is = Intersections { intersections: vec![i1, i2] };

    assert_eq!(is.hit(), None);
}

#[test]
fn hit_multiple() {
    let s  = ShapeNode::sphere();
    let i1 = Intersection { t: 5.0,  what: &s }; 
    let i2 = Intersection { t: 7.0,  what: &s };
    let i3 = Intersection { t: -3.0, what: &s };
    let i4 = Intersection { t: 2.0,  what: &s };
    let mut is = Intersections { intersections: vec![i1, i2, i3, i4] };

    assert_eq!(is.hit().unwrap(), i4);
}

#[test]
fn normal_on_sphere_x() {
    let s = ShapeNode::sphere();
    let n = normal_at(&s, Tuple4D::point(1.0, 0.0, 0.0));

    assert_eq!(n, Tuple4D::vector(1.0, 0.0, 0.0));
}

#[test]
fn normal_on_sphere_y() {
    let s = ShapeNode::sphere();
    let n = normal_at(&s, Tuple4D::point(0.0, 1.0, 0.0));

    assert_eq!(n, Tuple4D::vector(0.0, 1.0, 0.0));
}

#[test]
fn normal_on_sphere_z() {
    let s = ShapeNode::sphere();
    let n = normal_at(&s, Tuple4D::point(0.0, 0.0, 1.0));

    assert_eq!(n, Tuple4D::vector(0.0, 0.0, 1.0));
}

#[test]
fn normal_on_sphere_nonaxial() {
    let s = ShapeNode::sphere();
    let n = normal_at(&s,
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
    let mut s = ShapeNode::sphere();
    s.transform = Matrix4D::translation(0.0, 1.0, 0.0);
    let n = normal_at(&s, Tuple4D::point(0.0, 1.70711, -0.70711));
    
    assert_eq!(n, Tuple4D::vector(0.0, 0.70711, -0.70711));
}

#[test]
fn normal_on_sphere_transformed() {
    let mut s = ShapeNode::sphere();
    s.transform = Matrix4D::scaling(1.0, 0.5, 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 5.0);
    let n = normal_at(&s,
        Tuple4D::point(0.0, 2.0f64.sqrt() / 2.0, -(2.0f64.sqrt() / 2.0))
    );

    assert_eq!(n, Tuple4D::vector(0.0, 0.97014, -0.24254));
}

#[test]
fn precompute_intersection_state() {
    use crate::intersect::IntersectionComputation;

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let shape = ShapeNode::sphere();
    let i = Intersection { t: 4.0, what: &shape };

    let comps = IntersectionComputation::new(&r, &i, None);

    assert!(std::ptr::eq(comps.obj, i.what));
    assert_eq!(comps.t, i.t);
    assert_eq!(comps.point, Tuple4D::point(0.0, 0.0, -1.0));
    assert_eq!(comps.eyev, Tuple4D::vector(0.0, 0.0, -1.0));
    assert_eq!(comps.normalv, Tuple4D::vector(0.0, 0.0, -1.0));
}

#[test]
fn precompute_outside_intersection() {
    use crate::intersect::IntersectionComputation;

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let shape = ShapeNode::sphere();
    let i = Intersection { t: 4.0, what: &shape };

    let comps = IntersectionComputation::new(&r, &i, None);
    assert!(!comps.inside);
}

#[test]
fn precompute_inside_intersection() {
    use crate::intersect::IntersectionComputation;

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -0.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let shape = ShapeNode::sphere();
    let i = Intersection { t: 1.0, what: &shape };

    let comps = IntersectionComputation::new(&r, &i, None);

    assert!(comps.inside);
    assert!(std::ptr::eq(comps.obj, i.what));
    assert_eq!(comps.t, i.t);
    assert_eq!(comps.point, Tuple4D::point(0.0, 0.0, 1.0));
    assert_eq!(comps.eyev, Tuple4D::vector(0.0, 0.0, -1.0));
    assert_eq!(comps.normalv, Tuple4D::vector(0.0, 0.0, -1.0));
}

#[test]
fn hit_should_offset_point() {
    use crate::intersect::IntersectionComputation;

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let mut shape = ShapeNode::sphere();
    shape.transform = Matrix4D::translation(0.0, 0.0, 1.0);

    let i = Intersection { t: 5.0, what: &shape };
    let comps = IntersectionComputation::new(&r, &i, None);

    assert!(comps.over_point.z < -FEQ_EPSILON / 2.0);
    assert!(comps.point.z > comps.over_point.z);
}

#[test]
fn creating_a_shape_group() {
    let g = ShapeNode::group();
    let children = match g.ty {
        ShapeNodeType::Group(ref c) => c,
        _ => unreachable!(),
    };

    assert_eq!(g.transform, Matrix4D::identity());
    assert_eq!(children.len(), 0);
}

#[test]
fn shape_has_parent_attribute() {
    let s = ShapeNode::empty();
    let parent = s.parent.upgrade();

    assert!(parent.is_none());
}

#[test]
fn adding_a_child_to_a_shape_group() {
    let g = Rc::new(RefCell::new(ShapeNode::group()));
    let s = Rc::new(RefCell::new(ShapeNode::empty()));

    add_child_to_group(Rc::clone(&g), Rc::clone(&s));

    let gb = g.borrow();
    let children = match gb.ty {
        ShapeNodeType::Group(ref c) => c,
        _ => unreachable!(),
    };

    // Check that a child has been appended to group g.
    assert_eq!(children.len(), 1);

    // Make sure that the child has the same allocation as original ShapeNode.
    assert!(Rc::ptr_eq(&s, &children[0]));

    // Check that the child points to the parent Group.
    assert!(Weak::ptr_eq(&s.borrow().parent, &Rc::downgrade(&g)));
}

#[test]
fn intersecting_ray_with_empty_group() {
    let g = ShapeNode::group();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = g.local_intersect(&r);
    assert!(is.intersections.is_empty())
}

#[test]
fn intersecting_ray_with_nonempty_group() {
    let g = Rc::new(RefCell::new(ShapeNode::group()));
    let s1 = Rc::new(RefCell::new(ShapeNode::sphere()));
    let s2 = Rc::new(RefCell::new(ShapeNode::sphere()));
    s2.borrow_mut().transform = Matrix4D::translation(0.0, 0.0, -3.0);
    let s3 = Rc::new(RefCell::new(ShapeNode::sphere()));
    s3.borrow_mut().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    add_child_to_group(Rc::clone(&g), Rc::clone(&s1));
    add_child_to_group(Rc::clone(&g), Rc::clone(&s2));
    add_child_to_group(Rc::clone(&g), Rc::clone(&s3));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let gb = g.borrow();
    let is = gb.local_intersect(&r);
    assert_eq!(is.intersections.len(), 4);
    assert!(std::ptr::eq(is.intersections[0].what, Ref::leak(s2.borrow())));
    assert!(std::ptr::eq(is.intersections[1].what, Ref::leak(s2.borrow())));
    assert!(std::ptr::eq(is.intersections[2].what, Ref::leak(s1.borrow())));
    assert!(std::ptr::eq(is.intersections[3].what, Ref::leak(s1.borrow())));
}

#[test]
fn intersecting_a_transformed_group() {
    let g = Rc::new(RefCell::new(ShapeNode::group()));
    g.borrow_mut().transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let s1 = Rc::new(RefCell::new(ShapeNode::sphere()));
    s1.borrow_mut().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    add_child_to_group(Rc::clone(&g), Rc::clone(&s1));

    let r = Ray4D::new(
        Tuple4D::point(10.0, 0.0, -10.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let gb = g.borrow();
    let is = intersect(&gb, r);

    assert_eq!(is.intersections.len(), 2);
}

#[test]
fn converting_a_point_from_world_to_object_space() {
    let g1 = Rc::new(RefCell::new(ShapeNode::group()));
    g1.borrow_mut().transform = Matrix4D::rotation_y(std::f64::consts::PI / 2.);

    let g2 = Rc::new(RefCell::new(ShapeNode::group()));
    g2.borrow_mut().transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    // Group g2 is a child of group g1.
    add_child_to_group(Rc::clone(&g1), Rc::clone(&g2));

    let s = Rc::new(RefCell::new(ShapeNode::sphere()));
    s.borrow_mut().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    // Sphere s is a child of group g2.
    add_child_to_group(Rc::clone(&g2), Rc::clone(&s));

    let p = s.borrow().world_to_object(Tuple4D::point(-2.0, 0.0, -10.0));
    assert_eq!(p, Tuple4D::point(0.0, 0.0, -1.0));
}

#[test]
fn converting_a_normal_from_object_to_world_space() {
    let g1 = Rc::new(RefCell::new(ShapeNode::group()));
    g1.borrow_mut().transform = Matrix4D::rotation_y(std::f64::consts::PI / 2.);

    let g2 = Rc::new(RefCell::new(ShapeNode::group()));
    g2.borrow_mut().transform = Matrix4D::scaling(1.0, 2.0, 3.0);

    // Group g2 is a child of group g1.
    add_child_to_group(Rc::clone(&g1), Rc::clone(&g2));

    let s = Rc::new(RefCell::new(ShapeNode::sphere()));
    s.borrow_mut().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    // Sphere s is a child of group g2.
    add_child_to_group(Rc::clone(&g2), Rc::clone(&s));

    let normal = s.borrow().normal_to_world(Tuple4D::vector(
        3.0f64.sqrt() / 3.0, 3.0f64.sqrt() / 3.0, 3.0f64.sqrt() / 3.0
    ));

    assert_eq!(normal, Tuple4D::vector(0.2857, 0.4286, -0.8571));
}

#[test]
fn finding_the_normal_on_a_child_object() {
    let g1 = Rc::new(RefCell::new(ShapeNode::group()));
    g1.borrow_mut().transform = Matrix4D::rotation_y(std::f64::consts::PI / 2.);

    let g2 = Rc::new(RefCell::new(ShapeNode::group()));
    g2.borrow_mut().transform = Matrix4D::scaling(1.0, 2.0, 3.0);

    // Group g2 is a child of group g1.
    add_child_to_group(Rc::clone(&g1), Rc::clone(&g2));

    let s = Rc::new(RefCell::new(ShapeNode::sphere()));
    s.borrow_mut().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    // Sphere s is a child of group g2.
    add_child_to_group(Rc::clone(&g2), Rc::clone(&s));

    let normal = normal_at(
        &s.borrow(),
        Tuple4D::point(1.7321, 1.1547, -5.5774)
    );
    assert_eq!(normal, Tuple4D::vector(0.2857, 0.4286, -0.8571));
}
