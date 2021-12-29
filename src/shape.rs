use std::sync::{ Arc, Mutex, Weak };

use crate::consts::FEQ_EPSILON;
use crate::tuple::Tuple4D;
use crate::ray::Ray4D;
use crate::light::Material;
use crate::matrix::Matrix4D;
use crate::intersect::{ Intersection, Intersections, filter_intersections };
use crate::geometry::{ TriangleInfo, SmoothTriangleInfo, Bounds };

pub(crate) type ShapePtr = Arc<Mutex<ShapeNode>>;

#[derive(Debug)]
pub(crate) enum ShapeType {
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

    /// A triangle. See TriangleInfo for further explanation.
    // NOTE: As an aside, this is kind of a downside of using enums; it's not
    // easy to associate data with new variants.
    Triangle(TriangleInfo),

    /// A smooth triangle. See TriangleInfo for further explanation.
    SmoothTriangle(SmoothTriangleInfo),

    /// A group of shapes. Can include other groups of shapes.
    Group(Vec<ShapePtr>),

    /// A union of two shapes.
    Union(ShapePtr, ShapePtr),

    /// An intersection of two shapes.
    Intersection(ShapePtr, ShapePtr),

    /// A difference of two shapes.
    Difference(ShapePtr, ShapePtr),
}

/// Checks that two ShapeTypes are equal.
///
/// Value comparisons are made for each variant of `ShapeType` except for
/// `Group`, which compares the *pointers* of each group child.
///
/// See the `eq` implementation for more information.
impl PartialEq for ShapeType {
    fn eq(&self, other: &Self) -> bool {
        use ShapeType::*;
        match (self, other) {
            (Empty, Empty) => true,
            (Sphere, Sphere) => true,
            (Cube, Cube) => true,

            (Cone(lmin, lmax, lc), Cone(rmin, rmax, rc))
                => (lmin == rmin) && (lmax == rmax) && (lc == rc),

            (Cylinder(lmin, lmax, lc), Cylinder(rmin, rmax, rc))
                => (lmin == rmin) && (lmax == rmax) && (lc == rc),

            (Plane(ln), Plane(rn)) => ln == rn,
            (Triangle(lti), Triangle(rti)) => lti == rti,
            (SmoothTriangle(lsti), SmoothTriangle(rsti)) => lsti == rsti,

            (Group(ref lg), Group(ref rg)) => {
                if lg.len() != rg.len() {
                    false
                } else {
                    // Check that all children in the group are pointer-equal.
                    lg.iter().zip(rg.iter()).all(|(l, r)| {
                        Arc::ptr_eq(l, r)
                    })
                }
            },

            _ => false,
        }
    }
}

#[derive(Debug)]
pub(crate) struct ShapeNode {
    pub ty: ShapeType, 
    parent: Weak<Mutex<ShapeNode>>,

    pub transform: Matrix4D,
    pub material: Material,
    saved_bounds: Option<Bounds>,
}

impl Default for ShapeNode {
    fn default() -> ShapeNode {
        ShapeNode {
            ty: ShapeType::Empty,
            parent: Weak::new(),
            transform: Matrix4D::identity(),
            material: Default::default(),
            saved_bounds: None,
        }
    }
}

/// Checks that two ShapeNodes are equal.
///
/// Note that the `parent` field is not checked for equality; two equivalent
/// shapes can have different parents.
///
/// This may be fixed in a later commit.
impl PartialEq for ShapeNode {
    fn eq(&self, other: &Self) -> bool {
        self.ty == other.ty
            && self.transform == other.transform
            && self.material == other.material
    }
}

pub struct Shape {
    root: ShapePtr,
}

macro_rules! shape_constructor {
    ( $s:ident, $ty:expr ) => {
        pub fn $s() -> Shape {
            let shape = ShapeNode {
                ty: $ty,
                ..Default::default()
            };

            Shape {
                root: Arc::new(Mutex::new(shape))
            }
        }
    };
}

impl Shape {
    // Creates an empty shape which does nothing. Mostly for testing.
    shape_constructor!(empty, ShapeType::Empty);

    // Creates a unit sphere with identity transform and default material.
    shape_constructor!(sphere, ShapeType::Sphere);

    // Creates a plane with a normal pointing up along the Y axis.
    shape_constructor!(plane, ShapeType::Plane(Tuple4D::vector(0.0, 1.0, 0.0)));

    // Creates a unit cube with identity transform and default material.
    shape_constructor!(cube, ShapeType::Cube);

    // Creates an infinitely long cylinder with no end caps.
    shape_constructor!(cylinder,
        ShapeType::Cylinder(
            -1.0 * std::f64::INFINITY,
            std::f64::INFINITY,
            false
        )
    );

    // Creates a bounded cylinder without caps.
    shape_constructor!(bounded_cylinder,
        ShapeType::Cylinder(-1.0, 1.0, false));

    // Creates a bounded cylinder with caps.
    shape_constructor!(capped_cylinder,
        ShapeType::Cylinder(-1.0, 1.0, true));

    // Creates an bounded double-napped cone with no end caps.
    shape_constructor!(bounded_cone,
        ShapeType::Cone(-1.0, 1.0, false));

    // Creates a double-napped cone with end caps.
    shape_constructor!(capped_cone,
        ShapeType::Cone(-1.0, 1.0, true));

    /// Creates a triangle, defined by three points in space.
    pub fn triangle(p1: Tuple4D, p2: Tuple4D, p3: Tuple4D) -> ShapeNode {
        ShapeNode {
            ty: ShapeType::Triangle(TriangleInfo::new(p1, p2, p3)),
            ..Default::default()
        }
    }

    /// Creates a "smooth" triangle with normals at each vertex.
    pub fn smooth_triangle(p1: Tuple4D, p2: Tuple4D, p3: Tuple4D,
        n1: Tuple4D, n2: Tuple4D, n3: Tuple4D) -> ShapeNode {
        ShapeNode {
            ty: ShapeType::SmoothTriangle(
                SmoothTriangleInfo::new(p1, p2, p3, n1, n2, n3)
            ),
            ..Default::default()
        }
    }

    /// Creates a group, which holds a list of other shapes (possibly groups).
    pub fn group() -> ShapeNode {
        ShapeNode {
            ty: ShapeType::Group(Vec::new()),
            ..Default::default()
        }
    }

    pub fn csg_union(s1: ShapePtr, s2: ShapePtr) -> ShapePtr {
        let union_shape = Arc::new(Mutex::new(
            ShapeNode {
                // Temporarily set type to empty so that child shapes can have
                // some parent to refer to.
                ty: ShapeType::Empty,
                ..Default::default()
            }
        ));

        s1.lock().unwrap().parent = Arc::downgrade(&union_shape);
        s2.lock().unwrap().parent = Arc::downgrade(&union_shape);

        // Set the type to a union of the two shapes, then return the union.
        union_shape.lock().unwrap().ty = ShapeType::Union(s1, s2);
        union_shape
    }

    pub fn csg_intersection(s1: ShapePtr, s2: ShapePtr) -> ShapePtr {
        let intersection = Arc::new(Mutex::new(
            ShapeNode {
                // Temporarily set type to empty so that child shapes can have
                // some parent to refer to.
                ty: ShapeType::Empty,
                ..Default::default()
            }
        ));

        s1.lock().unwrap().parent = Arc::downgrade(&intersection);
        s2.lock().unwrap().parent = Arc::downgrade(&intersection);

        // Set the type to a intersection of the two shapes, then return it..
        intersection.lock().unwrap().ty = ShapeType::Intersection(s1, s2);
        intersection 
    }

    pub fn csg_difference(s1: ShapePtr, s2: ShapePtr) -> ShapePtr {
        let difference = Arc::new(Mutex::new(
            ShapeNode {
                // Temporarily set type to empty so that child shapes can have
                // some parent to refer to.
                ty: ShapeType::Empty,
                ..Default::default()
            }
        ));

        s1.lock().unwrap().parent = Arc::downgrade(&difference);
        s2.lock().unwrap().parent = Arc::downgrade(&difference);

        // Set the type to a difference of the two shapes, then return it..
        difference.lock().unwrap().ty = ShapeType::Difference(s1, s2);
        difference 
    }

    /// Gets the left operand of a CSG operator.
    ///
    /// Panics if `self` is not a CSG operator (types `Union`, `Intersection`
    /// and `Difference`).
    pub fn csg_left(&self) -> ShapePtr {
        Arc::clone(
            match self.ty {
                ShapeType::Union(ref l, _) => l,
                ShapeType::Intersection(ref l, _) => l,
                ShapeType::Difference(ref l, _) => l,
                _ => panic!("csg_left called on non-CSG type shape."),
            }
        )
    }

    /// Gets the right operand of a CSG operator.
    ///
    /// Panics if `self` is not a CSG operator (types `Union`, `Intersection`
    /// and `Difference`).
    pub fn csg_right(&self) -> ShapePtr {
        Arc::clone(
            match self.ty {
                ShapeType::Union(_, ref r) => r,
                ShapeType::Intersection(_, ref r) => r,
                ShapeType::Difference(_, ref r) => r,
                _ => panic!("csg_right called on non_CSG type shape."),
            }
        )
    }

    /// Deduces whether a shape includes another shape.
    ///
    /// If `self` is a `Group`, then this function will return true if any
    /// child of the `Group` satisfies `child.lock().unwrap().includes(other)`.
    ///
    /// If `self` is a CSG object (e.g. `Union`, `Intersection`), this function
    /// will return if either child of `self` includes `other`.
    ///
    /// If `self` is any other shape type, this function returns true if `self`
    /// is equal to `other`.
    pub fn includes(&self, other: &Self) -> bool {
        match self.ty {
            // Groups
            ShapeType::Group(ref children)
                => children.iter().any(|c| c.lock().unwrap().includes(other)),

            // CSG operations
            ShapeType::Union(ref left, ref right)
                => left.lock().unwrap().includes(other)
                    || right.lock().unwrap().includes(other),
            ShapeType::Intersection(ref left, ref right)
                => left.lock().unwrap().includes(other)
                    || right.lock().unwrap().includes(other),
            ShapeType::Difference(ref left, ref right)
                => left.lock().unwrap().includes(other)
                    || right.lock().unwrap().includes(other),

            // Other objects
            _ => self == other,
        }
    }

    /// Gets the bounds for different shape types.
    pub fn bounds(&self) -> Bounds {
        match self.ty {
            ShapeType::Empty => Bounds::empty(),
            ShapeType::Sphere => Bounds::new(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0),
            ShapeType::Cube => Bounds::new(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0),
            ShapeType::Cone(min, max, _)
                => Bounds::new(min, min, min, max, max, max),
            ShapeType::Cylinder(min, max, _)
                => Bounds::new(min, min, min, max, max, max),
            ShapeType::Triangle(ref ti) => {
                let min_x = ti.p1.x.min(ti.p2.x.min(ti.p3.x));
                let min_y = ti.p1.y.min(ti.p2.y.min(ti.p3.y));
                let min_z = ti.p1.z.min(ti.p2.z.min(ti.p3.z));
                let max_x = ti.p1.x.max(ti.p2.x.max(ti.p3.x));
                let max_y = ti.p1.y.max(ti.p2.y.max(ti.p3.y));
                let max_z = ti.p1.z.max(ti.p2.z.max(ti.p3.z));

                Bounds::new(min_x, min_y, min_z, max_x, max_y, max_z)
            },
            ShapeType::SmoothTriangle(ref sti) => {
                let p1 = sti.triangle_info.p1;
                let p2 = sti.triangle_info.p2;
                let p3 = sti.triangle_info.p3;

                let min_x = p1.x.min(p2.x.min(p3.x));
                let min_y = p1.y.min(p2.y.min(p3.y));
                let min_z = p1.z.min(p2.z.min(p3.z));
                let max_x = p1.x.max(p2.x.max(p3.x));
                let max_y = p1.y.max(p2.y.max(p3.y));
                let max_z = p1.z.max(p2.z.max(p3.z));

                Bounds::new(min_x, min_y, min_z, max_x, max_y, max_z)
            },
            ShapeType::Group(ref children) => {
                // Collect the bounds of each child.
                let child_bounds: Vec<Bounds> = children.iter().map(|child|
                    // Transform each child bounds into group space.
                    child.lock().unwrap().bounds().transform(&self.transform)
                ).collect();

                let mut min_x = std::f64::INFINITY;
                let mut min_y = std::f64::INFINITY;
                let mut min_z = std::f64::INFINITY;
                let mut max_x = -1.0 * std::f64::INFINITY;
                let mut max_y = -1.0 * std::f64::INFINITY;
                let mut max_z = -1.0 * std::f64::INFINITY;

                // Using the child bounds, find the minimum and maximum bounds
                // which will encapsulate all child bounds.
                for bounds in child_bounds {
                    min_x = min_x.min(bounds.minimum.x);
                    min_y = min_y.min(bounds.minimum.y);
                    min_z = min_z.min(bounds.minimum.z);
                    max_x = max_x.max(bounds.maximum.x);
                    max_y = max_y.max(bounds.maximum.y);
                    max_z = max_z.max(bounds.maximum.z);
                }

                Bounds::new(min_x, min_y, min_z, max_x, max_y, max_z)
            },

            // TODO: How to handle planes? Should we assume that normal always
            // points up?
            _ => Default::default(),
        }
    }

    /// Returns a reference to a list of child `ShapeNode`s if this is a group.
    pub fn children(&self) -> Option<&Vec<ShapePtr>> {
        if let ShapeType::Group(ref children) = self.ty {
            Some(children)
        } else {
            None
        }
    }

    /// Returns a reference to the `TriangleInfo` if this is a triangle.
    pub fn triangle_info(&self) -> Option<&TriangleInfo> {
        if let ShapeType::Triangle(ref info) = self.ty {
            Some(info)
        } else {
            None
        }
    }

    pub fn smooth_triangle_info(&self) -> Option<&SmoothTriangleInfo> {
        if let ShapeType::SmoothTriangle(ref smooth_info) = self.ty {
            Some(smooth_info)
        } else {
            None
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
            point = parent.lock().unwrap().world_to_object(point);
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
            normal = parent.lock().unwrap().normal_to_world(normal);
        }

        normal
    }

    pub fn local_intersect(&self, ray: &Ray4D) -> Intersections {
        match self.ty {
            ShapeType::Empty => Intersections::new(),
            ShapeType::Sphere => self.intersect_sphere(ray),
            ShapeType::Plane(_) => self.intersect_plane(ray),
            ShapeType::Cube => self.intersect_cube(ray),
            ShapeType::Cylinder(_, _, _) => self.intersect_cylinder(ray),
            ShapeType::Cone(_, _, _) => self.intersect_cone(ray),
            ShapeType::Triangle(_) => self.intersect_triangle(ray),
            ShapeType::SmoothTriangle(_) => self.intersect_smooth_triangle(ray),
            ShapeType::Group(_) => self.intersect_group(ray),
            ShapeType::Union(_,_)
                | ShapeType::Intersection(_,_)
                | ShapeType::Difference(_,_) => self.intersect_csg(ray),
        }
    }

    pub fn local_normal_at(&self, at: &Tuple4D, hit: &Intersection) -> Tuple4D {
        match self.ty {
            ShapeType::Empty => Tuple4D { w: 0.0, ..*at },
            ShapeType::Sphere => self.normal_at_sphere(at),
            ShapeType::Plane(_) => self.normal_at_plane(at),
            ShapeType::Cube => self.normal_at_cube(at),
            ShapeType::Cylinder(_, _, _) => self.normal_at_cylinder(at),
            ShapeType::Cone(_, _, _) => self.normal_at_cone(at),
            ShapeType::Triangle(_) => self.normal_at_triangle(at),
            ShapeType::SmoothTriangle(_)
                => self.normal_at_smooth_triangle(at, hit),
            ShapeType::Group(_) => panic!(
                "Local normal calculations should never occur on groups."
            ),

            // CSG operations always defer intersections to their children
            // objects, so no normals will ever be calculated on CSGs.
            ShapeType::Union(_,_)
                | ShapeType::Intersection(_,_)
                | ShapeType::Difference(_,_) => panic!(
                "Local normal calculations should never occur on CSG ops."
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
            ShapeType::Sphere => (),
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

        let i1 = Intersection::new(t1, self);
        let i2 = Intersection::new(t2, self);
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
        let i = Intersection::new(t, self);

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
                Intersection::new(tmin, self),
                Intersection::new(tmax, self) 
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
                Intersection::new(t0, self)
            );
        }

        // If the t1 intersection is within the cylinder's bounds, add it.
        let y1 = ray.origin.y + t1 * ray.direction.y;
        if minimum < y1 && y1 < maximum {
            is.intersections.push(
                Intersection::new(t1, self)
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
                is.intersections.push(Intersection::new(t, self));
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
                Intersection::new(t0, self)
            );
        }

        // If the t1 intersection is within the cone's bounds, add it.
        let y1 = ray.origin.y + t1 * ray.direction.y;
        if minimum < y1 && y1 < maximum {
            is.intersections.push(
                Intersection::new(t1, self)
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

    fn intersect_group(&self, ray: &Ray4D) -> Intersections {
        let children = match self.ty {
            ShapeType::Group(ref c) => c,
            _ => unreachable!(),
        };

        // If there are no children, return an empty list of intersections.
        if children.is_empty() {
            return Intersections::new()
        }

        // Get the bounds to see if anything in this group will be intersected.
        let bounds = self.bounds();

        let (xtmin, xtmax) = Bounds::check_axis(
            bounds.minimum.x, bounds.maximum.x, ray.origin.x, ray.direction.x
        );
        let (ytmin, ytmax) = Bounds::check_axis(
            bounds.minimum.y, bounds.maximum.y, ray.origin.y, ray.direction.y
        );
        let (ztmin, ztmax) = Bounds::check_axis(
            bounds.minimum.z, bounds.maximum.z, ray.origin.z, ray.direction.z
        );

        // Calculate the minimum t and maximum t along the bounds box.
        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);

        // If the ray doesn't hit the bounds box, return no intersections for
        // this group. (No children can be hit).
        if tmin > tmax {
            return Intersections::new()
        }

        // Otherwise, for each child, collect its intersections and aggregate.
        let mut all_intersections = Vec::new();
        for child in children.iter() {
            // Use regular intersect, otherwise transforms won't be applied
            let is = intersect(&child.lock().unwrap(), *ray);
            all_intersections.push(is)
        }

        Intersections::aggregate(all_intersections)
    }

    fn intersect_triangle(&self, ray: &Ray4D) -> Intersections {
        let triangle_info = match self.ty {
            ShapeType::Triangle(ref ti) => ti,
            _ => unreachable!(),
        };

        let dir_cross_e2 = ray.direction.cross(&triangle_info.e2);
        let determinant = triangle_info.e1.dot(&dir_cross_e2);

        // If the ray is parallel to the triangle, return no intersections.
        if determinant.abs() < FEQ_EPSILON {
            return Intersections::new();
        }

        let f = 1.0 / determinant;
        let p1_to_origin = ray.origin - triangle_info.p1;
        let u = f * p1_to_origin.dot(&dir_cross_e2);
        if u < 0.0 {
            return Intersections::new();
        }

        let origin_cross_e1 = p1_to_origin.cross(&triangle_info.e1);
        let v = f * ray.direction.dot(&origin_cross_e1);
        if v < 0.0 || u + v > 1.0 {
            return Intersections::new();
        }

        // If an intersection has occurred, return it.
        // Note that triangles can only be intersected on one side.
        let t = f * triangle_info.e2.dot(&origin_cross_e1);
        Intersections {
            intersections: vec![
                Intersection::new(t, &self)
            ],
        }
    }

    fn normal_at_triangle(&self, _at: &Tuple4D) -> Tuple4D {
        let triangle_info = match self.ty {
            ShapeType::Triangle(ref ti) => ti,
            _ => unreachable!(),
        };

        triangle_info.normal
    }

    fn intersect_smooth_triangle(&self, ray: &Ray4D) -> Intersections {
        // TODO reduce code duplication with `intersect_triangle`

        let smooth_triangle_info = match self.ty {
            ShapeType::SmoothTriangle(ref sti) => sti,
            _ => unreachable!(),
        };

        let dir_cross_e2 = ray.direction.cross(
            &smooth_triangle_info.triangle_info.e2
        );
        let determinant =
            smooth_triangle_info.triangle_info.e1.dot(&dir_cross_e2);

        // If the ray is parallel to the triangle, return no intersections.
        if determinant.abs() < FEQ_EPSILON {
            return Intersections::new();
        }

        let f = 1.0 / determinant;
        let p1_to_origin = ray.origin - smooth_triangle_info.triangle_info.p1;
        let u = f * p1_to_origin.dot(&dir_cross_e2);
        if u < 0.0 {
            return Intersections::new();
        }

        let origin_cross_e1 = p1_to_origin.cross(
            &smooth_triangle_info.triangle_info.e1
        );
        let v = f * ray.direction.dot(&origin_cross_e1);
        if v < 0.0 || u + v > 1.0 {
            return Intersections::new();
        }

        // If an intersection has occurred, return it.
        // Note that triangles can only be intersected on one side.
        let t = f * smooth_triangle_info.triangle_info.e2.dot(&origin_cross_e1);
        Intersections {
            intersections: vec![
                Intersection::new_uv(t, &self, u, v)
            ],
        }
    }

    fn normal_at_smooth_triangle(&self, _at: &Tuple4D, hit: &Intersection)
        -> Tuple4D {
        let smooth_triangle_info = match self.ty {
            ShapeType::SmoothTriangle(ref sti) => sti,
            _ => unreachable!(),
        };

        if let Some((u,v)) = hit.uv {
            smooth_triangle_info.n2 * u
                + smooth_triangle_info.n3 * v
                + smooth_triangle_info.n1 * (1.0 - u - v)
        } else {
            smooth_triangle_info.triangle_info.normal
        }
    }

    /// Intersects a CSG operator.
    ///
    /// Each CSG operator gives two operands, of which each holds a pointer to
    /// some child shapes. These shapes are coalesced as expected through each
    /// CSG operation.
    fn intersect_csg(&self, ray: &Ray4D) -> Intersections {
        let (left_csg, right_csg) = match self.ty {
            ShapeType::Union(ref lc, ref rc) => (lc, rc),
            ShapeType::Intersection(ref lc, ref rc) => (lc, rc),
            ShapeType::Difference(ref lc, ref rc) => (lc, rc),
            _ => unreachable!(),
        };

        // Get the intersection of the left/right children
        let left_is = intersect(&left_csg.lock().unwrap(), *ray);
        let right_is = intersect(&right_csg.lock().unwrap(), *ray);

        // Function aggregate sorts all intersections provided.
        let all_is = Intersections::aggregate(vec![left_is, right_is]);

        filter_intersections(self, &all_is) 
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

    fn intersect_cylinder_caps(&self, ray: &Ray4D, is: &mut Intersections) {
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
            is.intersections.push(Intersection::new(tl, self));
        }

        // Check for an intersection with the upper end cap.
        let tu = (maximum - ray.origin.y) / ray.direction.y;
        if Self::check_cylinder_cap(ray, tu) {
            is.intersections.push(Intersection::new(tu, self));
        }
    }

    /// Checks to see if the intersection `t` is within a radius of 1 (the
    /// assumed radius of a cylinder) from the Y axis.
    fn check_cylinder_cap(ray: &Ray4D, t: f64) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;

        (x.powi(2) + z.powi(2)) <= 1.0
    }

    fn intersect_cone_caps(&self, ray: &Ray4D, is: &mut Intersections) {
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
            is.intersections.push(Intersection::new(tl, self));
        }

        // Check for an intersection with the upper end cap.
        let tu = (maximum - ray.origin.y) / ray.direction.y;
        if Self::check_cone_cap(ray, tu, maximum) {
            is.intersections.push(Intersection::new(tu, self));
        }
    }

    fn check_cone_cap(ray: &Ray4D, t: f64, y: f64) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;

        x.powi(2) + z.powi(2) <= y.powi(2)
    }
}

pub fn add_child_to_group(group: ShapePtr, child: ShapePtr) {
    let mut gb = group.lock().unwrap();
    let children = match gb.ty {
        ShapeType::Group(ref mut c) => c,
        _ => panic!("Cannot add child to non-group shape."),
    };

    child.lock().unwrap().parent = Arc::downgrade(&group);
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

pub fn normal_at(s: &ShapeNode, world_point: Tuple4D, hit: &Intersection)
    -> Tuple4D {
    let local_point = s.world_to_object(world_point);
    let local_normal = s.local_normal_at(&local_point, hit);
    s.normal_to_world(local_normal)
}

#[test]
fn compute_normal_on_translated_sphere() {
    let mut s = ShapeNode::sphere();
    s.transform = Matrix4D::translation(0.0, 1.0, 0.0);

    let p = Tuple4D::point(0.0, 1.70711, -0.70711);
    let n = normal_at(&s, p, &Intersection::new(0.0, &s));

    assert_eq!(n, Tuple4D::vector(0.0, 0.70711, -0.70711));
}

#[test]
fn compute_normal_on_transformed_sphere() {
    let mut s = ShapeNode::sphere();
    s.transform = Matrix4D::scaling(1.0, 0.5, 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 5.0);
    
    let p = Tuple4D::point(0.0, 2.0f64.sqrt() / 2.0, -(2.0f64.sqrt()) / 2.0);
    let n = normal_at(&s, p, &Intersection::new(0.0, &s));

    assert_eq!(n, Tuple4D::vector(0.0, 0.97014, -0.24254));
}

#[test]
fn normal_on_plane() {
    let p = ShapeNode::plane();

    let n1 = p.local_normal_at(
        &Tuple4D::point(0.0, 0.0, 0.0),
        &Intersection::new(0.0, &p)
    );
    let n2 = p.local_normal_at(
        &Tuple4D::point(10.0, 0.0, -10.0),
        &Intersection::new(0.0, &p)
    );
    let n3 = p.local_normal_at(
        &Tuple4D::point(-5.0, 0.0, 150.0),
        &Intersection::new(0.0, &p)
    );

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
    let i1 = Intersection::new(1.0, &s); 
    let i2 = Intersection::new(2.0, &s);
    let mut is = Intersections { intersections: vec![i1, i2] };

    assert_eq!(is.hit().unwrap(), i1);
}

#[test]
fn hit_with_some_negative() {
    let s  = ShapeNode::sphere();
    let i1 = Intersection::new(-1.0, &s); 
    let i2 = Intersection::new( 1.0, &s);
    let mut is = Intersections { intersections: vec![i1, i2] };

    assert_eq!(is.hit().unwrap(), i2);
}

#[test]
fn hit_with_all_negative() {
    let s  = ShapeNode::sphere();
    let i1 = Intersection::new(-2.0, &s);
    let i2 = Intersection::new(-1.0, &s);
    let mut is = Intersections { intersections: vec![i1, i2] };

    assert_eq!(is.hit(), None);
}

#[test]
fn hit_multiple() {
    let s  = ShapeNode::sphere();
    let i1 = Intersection::new(5.0,  &s); 
    let i2 = Intersection::new(7.0,  &s);
    let i3 = Intersection::new(-3.0, &s);
    let i4 = Intersection::new(2.0,  &s);
    let mut is = Intersections { intersections: vec![i1, i2, i3, i4] };

    assert_eq!(is.hit().unwrap(), i4);
}

#[test]
fn normal_on_sphere_x() {
    let s = ShapeNode::sphere();
    let n = normal_at(
        &s, Tuple4D::point(1.0, 0.0, 0.0), &Intersection::new(0.0, &s)
    );

    assert_eq!(n, Tuple4D::vector(1.0, 0.0, 0.0));
}

#[test]
fn normal_on_sphere_y() {
    let s = ShapeNode::sphere();
    let n = normal_at(
        &s, Tuple4D::point(0.0, 1.0, 0.0), &Intersection::new(0.0, &s)
    );

    assert_eq!(n, Tuple4D::vector(0.0, 1.0, 0.0));
}

#[test]
fn normal_on_sphere_z() {
    let s = ShapeNode::sphere();
    let n = normal_at(
        &s, Tuple4D::point(0.0, 0.0, 1.0), &Intersection::new(0.0, &s)
    );

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
        ),
        &Intersection::new(0.0, &s)
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
    let n = normal_at(
        &s, Tuple4D::point(0.0, 1.70711, -0.70711), &Intersection::new(0.0, &s)
    );
    
    assert_eq!(n, Tuple4D::vector(0.0, 0.70711, -0.70711));
}

#[test]
fn normal_on_sphere_transformed() {
    let mut s = ShapeNode::sphere();
    s.transform = Matrix4D::scaling(1.0, 0.5, 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 5.0);
    let n = normal_at(&s,
        Tuple4D::point(0.0, 2.0f64.sqrt() / 2.0, -(2.0f64.sqrt() / 2.0)),
        &Intersection::new(0.0, &s)
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
    let i = Intersection::new(4.0, &shape);

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
    let i = Intersection::new(4.0, &shape);

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
    let i = Intersection::new(1.0, &shape);

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

    let i = Intersection::new(5.0, &shape);
    let comps = IntersectionComputation::new(&r, &i, None);

    assert!(comps.over_point.z < -FEQ_EPSILON / 2.0);
    assert!(comps.point.z > comps.over_point.z);
}

#[test]
fn creating_a_shape_group() {
    let g = ShapeNode::group();
    let children = match g.ty {
        ShapeType::Group(ref c) => c,
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
    let g = Arc::new(Mutex::new(ShapeNode::group()));
    let s = Arc::new(Mutex::new(ShapeNode::empty()));

    add_child_to_group(Arc::clone(&g), Arc::clone(&s));

    let gb = g.lock().unwrap();
    let children = match gb.ty {
        ShapeType::Group(ref c) => c,
        _ => unreachable!(),
    };

    // Check that a child has been appended to group g.
    assert_eq!(children.len(), 1);

    // Make sure that the child has the same allocation as original ShapeNode.
    assert!(Rc::ptr_eq(&s, &children[0]));

    // Check that the child points to the parent Group.
    assert!(Weak::ptr_eq(&s.lock().unwrap().parent, &Arc::downgrade(&g)));
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
    let g = Arc::new(Mutex::new(ShapeNode::group()));
    let s1 = Arc::new(Mutex::new(ShapeNode::sphere()));
    let s2 = Arc::new(Mutex::new(ShapeNode::sphere()));
    s2.lock().unwrap().transform = Matrix4D::translation(0.0, 0.0, -3.0);
    let s3 = Arc::new(Mutex::new(ShapeNode::sphere()));
    s3.lock().unwrap().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    add_child_to_group(Arc::clone(&g), Arc::clone(&s1));
    add_child_to_group(Arc::clone(&g), Arc::clone(&s2));
    add_child_to_group(Arc::clone(&g), Arc::clone(&s3));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let gb = g.lock().unwrap();
    let is = gb.local_intersect(&r);
    assert_eq!(is.intersections.len(), 4);
    assert!(std::ptr::eq(is.intersections[0].what, Ref::leak(s2.lock().unwrap())));
    assert!(std::ptr::eq(is.intersections[1].what, Ref::leak(s2.lock().unwrap())));
    assert!(std::ptr::eq(is.intersections[2].what, Ref::leak(s1.lock().unwrap())));
    assert!(std::ptr::eq(is.intersections[3].what, Ref::leak(s1.lock().unwrap())));
}

#[test]
fn intersecting_a_transformed_group() {
    let g = Arc::new(Mutex::new(ShapeNode::group()));
    g.lock().unwrap().transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    let s1 = Arc::new(Mutex::new(ShapeNode::sphere()));
    s1.lock().unwrap().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    add_child_to_group(Arc::clone(&g), Arc::clone(&s1));

    let r = Ray4D::new(
        Tuple4D::point(10.0, 0.0, -10.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let gb = g.lock().unwrap();
    let is = intersect(&gb, r);

    assert_eq!(is.intersections.len(), 2);
}

#[test]
fn converting_a_point_from_world_to_object_space() {
    let g1 = Arc::new(Mutex::new(ShapeNode::group()));
    g1.lock().unwrap().transform = Matrix4D::rotation_y(std::f64::consts::PI / 2.);

    let g2 = Arc::new(Mutex::new(ShapeNode::group()));
    g2.lock().unwrap().transform = Matrix4D::scaling(2.0, 2.0, 2.0);

    // Group g2 is a child of group g1.
    add_child_to_group(Arc::clone(&g1), Arc::clone(&g2));

    let s = Arc::new(Mutex::new(ShapeNode::sphere()));
    s.lock().unwrap().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    // Sphere s is a child of group g2.
    add_child_to_group(Arc::clone(&g2), Arc::clone(&s));

    let p = s.lock().unwrap().world_to_object(Tuple4D::point(-2.0, 0.0, -10.0));
    assert_eq!(p, Tuple4D::point(0.0, 0.0, -1.0));
}

#[test]
fn converting_a_normal_from_object_to_world_space() {
    let g1 = Arc::new(Mutex::new(ShapeNode::group()));
    g1.lock().unwrap().transform = Matrix4D::rotation_y(std::f64::consts::PI / 2.);

    let g2 = Arc::new(Mutex::new(ShapeNode::group()));
    g2.lock().unwrap().transform = Matrix4D::scaling(1.0, 2.0, 3.0);

    // Group g2 is a child of group g1.
    add_child_to_group(Arc::clone(&g1), Arc::clone(&g2));

    let s = Arc::new(Mutex::new(ShapeNode::sphere()));
    s.lock().unwrap().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    // Sphere s is a child of group g2.
    add_child_to_group(Arc::clone(&g2), Arc::clone(&s));

    let normal = s.lock().unwrap().normal_to_world(Tuple4D::vector(
        3.0f64.sqrt() / 3.0, 3.0f64.sqrt() / 3.0, 3.0f64.sqrt() / 3.0
    ));

    assert_eq!(normal, Tuple4D::vector(0.2857, 0.4286, -0.8571));
}

#[test]
fn finding_the_normal_on_a_child_object() {
    let g1 = Arc::new(Mutex::new(ShapeNode::group()));
    g1.lock().unwrap().transform = Matrix4D::rotation_y(std::f64::consts::PI / 2.);

    let g2 = Arc::new(Mutex::new(ShapeNode::group()));
    g2.lock().unwrap().transform = Matrix4D::scaling(1.0, 2.0, 3.0);

    // Group g2 is a child of group g1.
    add_child_to_group(Arc::clone(&g1), Arc::clone(&g2));

    let s = Arc::new(Mutex::new(ShapeNode::sphere()));
    s.lock().unwrap().transform = Matrix4D::translation(5.0, 0.0, 0.0);

    // Sphere s is a child of group g2.
    add_child_to_group(Arc::clone(&g2), Arc::clone(&s));

    let normal = normal_at(
        &s.lock().unwrap(),
        Tuple4D::point(1.7321, 1.1547, -5.5774),
        &Intersection::new(0.0, &ShapeNode::sphere())
    );
    assert_eq!(normal, Tuple4D::vector(0.2857, 0.4286, -0.8571));
}

#[test]
fn constructing_a_triangle() {
    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let t = ShapeNode::triangle(p1, p2, p3);

    if let ShapeType::Triangle(ti) = t.ty {
        assert_eq!(ti.p1, p1);
        assert_eq!(ti.p2, p2);
        assert_eq!(ti.p3, p3);
        assert_eq!(ti.e1, Tuple4D::vector(-1.0, -1.0, 0.0));
        assert_eq!(ti.e2, Tuple4D::vector(1.0, -1.0, 0.0));
        assert_eq!(ti.normal, Tuple4D::vector(0.0, 0.0, -1.0));
    } else {
        unreachable!();
    }
}

#[test]
fn finding_the_normal_on_a_triangle() {
    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let t = ShapeNode::triangle(p1, p2, p3);

    let i = Intersection::new(0.0, &t);
    if let ShapeType::Triangle(ref ti) = t.ty {
        assert_eq!(ti.normal,
            normal_at(&t, Tuple4D::point(0.0, 0.5, 0.0), &i));
        assert_eq!(ti.normal,
            normal_at(&t, Tuple4D::point(-0.5, 0.75, 0.0), &i));
        assert_eq!(ti.normal,
            normal_at(&t, Tuple4D::point(0.5, 0.25, 0.0), &i));
    } else {
        unreachable!();
    }
}

#[test]
fn intersecting_a_ray_parallel_to_a_triangle() {
    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let t = ShapeNode::triangle(p1, p2, p3);

    let r = Ray4D::new(
        Tuple4D::point(0.0, -1.0, -2.0),
        Tuple4D::vector(0.0, 1.0, 0.0)
    );

    let is = t.intersect_triangle(&r);
    assert!(is.intersections.is_empty())
}

#[test]
fn a_ray_misses_the_p1_p3_edge() {
    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let t = ShapeNode::triangle(p1, p2, p3);   

    let r = Ray4D::new(
        Tuple4D::point(1.0, 1.0, -2.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = t.intersect_triangle(&r);
    assert!(is.intersections.is_empty())
}

#[test]
fn a_ray_misses_the_p1_p2_edge() {
    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let t = ShapeNode::triangle(p1, p2, p3);   

    let r = Ray4D::new(
        Tuple4D::point(-1.0, 1.0, -2.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = t.intersect_triangle(&r);
    assert!(is.intersections.is_empty())
}

#[test]
fn a_ray_misses_the_p2_p3_edge() {
    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let t = ShapeNode::triangle(p1, p2, p3);   

    let r = Ray4D::new(
        Tuple4D::point(0.0, -1.0, -2.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = t.intersect_triangle(&r);
    assert!(is.intersections.is_empty())
}

#[test]
fn a_ray_strikes_a_triangle() {
    use crate::feq;

    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let t = ShapeNode::triangle(p1, p2, p3);   

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.5, -2.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = t.intersect_triangle(&r);
    assert_eq!(is.intersections.len(), 1);
    assert!(feq(is.intersections[0].t, 2.0));
}

#[test]
fn constructing_a_smooth_triangle() {
    let p1 = Tuple4D::point(0.0, 1.0, 0.0);
    let p2 = Tuple4D::point(-1.0, 0.0, 0.0);
    let p3 = Tuple4D::point(1.0, 0.0, 0.0);
    let n1 = Tuple4D::vector(1.0, 0.0, 0.0);
    let n2 = Tuple4D::vector(-1.0, 0.0, 0.0);
    let n3 = Tuple4D::vector(1.0, 0.0, 0.0);

    let s = ShapeNode::smooth_triangle(p1, p2, p3, n1, n2, n3);
    let smooth_triangle_info = match s.ty {
        ShapeType::SmoothTriangle(sti) => sti,
        _ => unreachable!(),
    };

    assert_eq!(smooth_triangle_info.triangle_info.p1, p1);
    assert_eq!(smooth_triangle_info.triangle_info.p2, p2);
    assert_eq!(smooth_triangle_info.triangle_info.p3, p3);
    assert_eq!(smooth_triangle_info.n1, n1);
    assert_eq!(smooth_triangle_info.n2, n2);
    assert_eq!(smooth_triangle_info.n3, n3);
}

#[test]
fn an_intersection_can_encapsulate_u_and_v() {
    let s = ShapeNode::triangle(
        Tuple4D::point(0.0, 1.0, 0.0),
        Tuple4D::point(-1.0, 0.0, 0.0),
        Tuple4D::point(1.0, 0.0, 0.0)
    );

    let i = Intersection::new_uv(3.5, &s, 0.2, 0.4);
    
    if let Some((u, v)) = i.uv {
        assert_eq!(u, 0.2);
        assert_eq!(v, 0.4);
    } else {
        unreachable!();
    }
}

#[test]
fn an_intersection_with_a_smooth_triangle_stores_uv() {
    use crate::feq;

    let s = ShapeNode::smooth_triangle(
        Tuple4D::point(0.0, 1.0, 0.0),
        Tuple4D::point(-1.0, 0.0, 0.0),
        Tuple4D::point(1.0, 0.0, 0.0),
        Tuple4D::vector(1.0, 0.0, 0.0),
        Tuple4D::vector(-1.0, 0.0, 0.0),
        Tuple4D::vector(1.0, 0.0, 0.0)
    );

    let r = Ray4D::new(
        Tuple4D::point(-0.2, 0.3, -2.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = s.local_intersect(&r);
    if let Some((u, v)) = is.intersections[0].uv {
        assert!(feq(u, 0.45));
        assert!(feq(v, 0.25));
    } else {
        unreachable!();
    }
}

#[test]
fn a_smooth_triangle_uses_uv_to_interpolate_the_normal() {
    let s = ShapeNode::smooth_triangle(
        Tuple4D::point(0.0, 1.0, 0.0),
        Tuple4D::point(-1.0, 0.0, 0.0),
        Tuple4D::point(1.0, 0.0, 0.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
        Tuple4D::vector(-1.0, 0.0, 0.0),
        Tuple4D::vector(1.0, 0.0, 0.0)
    );

    let i = Intersection::new_uv(1.0, &s, 0.45, 0.25);
    let n = normal_at(&s, Tuple4D::point(0.0, 0.0, 0.0), &i);

    assert_eq!(n, Tuple4D::vector(-0.5547, 0.83205, 0.0));
}

#[test]
fn csg_is_created_with_an_operation_and_two_shapes() {
    let s1 = Arc::new(Mutex::new(ShapeNode::sphere()));
    let s2 = Arc::new(Mutex::new(ShapeNode::cube()));

    let c = ShapeNode::csg_union(Arc::clone(&s1), Arc::clone(&s2));
    let c_ref = c.lock().unwrap();
    if let ShapeType::Union(ref s1_ptr, ref s2_ptr) = c_ref.ty {
        assert!(Rc::ptr_eq(s1_ptr, &s1));
        assert!(Rc::ptr_eq(s2_ptr, &s2));
        assert!(Weak::ptr_eq(&Arc::downgrade(&c), &s1.lock().unwrap().parent));
        assert!(Weak::ptr_eq(&Arc::downgrade(&c), &s2.lock().unwrap().parent));
    }
}

#[test]
fn a_ray_misses_a_csg_object() {
    let s1 = Arc::new(Mutex::new(ShapeNode::sphere()));
    let s2 = Arc::new(Mutex::new(ShapeNode::cube()));
    let c = ShapeNode::csg_union(s1, s2);

    let r = Ray4D::new(
        Tuple4D::point(0.0, 2.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let c_ref = c.lock().unwrap();
    let is = c_ref.local_intersect(&r);
    assert!(is.intersections.is_empty());
}

#[test]
fn a_ray_hits_a_csg_object() {
    use crate::feq;

    let s1 = Arc::new(Mutex::new(ShapeNode::sphere()));
    let s2 = Arc::new(Mutex::new(ShapeNode::sphere()));
    s2.lock().unwrap().transform = Matrix4D::translation(0.0, 0.0, 0.5);

    let c = ShapeNode::csg_union(Arc::clone(&s1), Arc::clone(&s2));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let c_ref = c.lock().unwrap();
    let is = c_ref.local_intersect(&r);

    assert_eq!(is.intersections.len(), 2);
    assert!(feq(is.intersections[0].t, 4.0));
    assert!(std::ptr::eq(is.intersections[0].what, Ref::leak(s1.lock().unwrap())));
    assert!(feq(is.intersections[1].t, 6.5));
    assert!(std::ptr::eq(is.intersections[1].what, Ref::leak(s2.lock().unwrap())));
}
