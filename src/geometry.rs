use std::rc::Rc;
use std::cell::RefCell;

use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::ray::Ray4D;
use crate::light::Material;

#[derive(Clone, Debug)]
pub struct Intersection {
    pub t: f64,
    pub what: Rc<RefCell<dyn IntersectableDebug>>,
}

impl PartialEq for Intersection {
    fn eq(&self, other: &Intersection) -> bool {
        self.t == other.t && Rc::ptr_eq(&self.what, &other.what)
    }
}

pub struct Intersections {
    pub intersections: Vec<Intersection>,
}

impl Intersections {
    pub fn new() -> Intersections {
        Intersections { intersections: Vec::new() }
    }

    pub fn hit<'a>(&'a mut self) -> Option<&'a Intersection> {
        self.intersections.retain(|i| i.t.is_finite());
        self.intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());

        for i in self.intersections.iter() {
            if i.t >= 0.0 {
                return Some(i);
            }
        }

        None
    }

    pub fn sort(&mut self) {
        self.intersections.sort_by(|a, b|
            a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal)
        );
    }
}

pub trait IntersectableDebug : Intersectable + std::fmt::Debug { }

pub trait Intersectable {
    fn intersect(&self, ray: Ray4D) -> Intersections;
    fn normal(&self, at: Tuple4D) -> Tuple4D;
    fn material(&self) -> &Material;
    fn material_mut(&mut self) -> &mut Material;
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Sphere {
    pub pos: Tuple4D,
    pub radius: f64,

    pub transform: Matrix4D,
    pub material: Material,
}

impl Sphere {
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

    pub fn unit() -> Sphere {
        Sphere {
            pos: Tuple4D::point(0.0, 0.0, 0.0),
            radius: 1.0,
            transform: Matrix4D::identity(),
            material: Default::default()
        }
    }
}

impl Intersectable for Sphere {
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

    fn material(&self) -> &Material {
        &self.material
    }

    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }
}

impl IntersectableDebug for Sphere { }

#[derive(Clone, Debug)]
pub struct IntersectionComputation {
    pub t: f64,
    pub obj: Rc<RefCell<dyn IntersectableDebug>>,

    pub point: Tuple4D,
    pub eyev: Tuple4D,
    pub normalv: Tuple4D,

    pub inside: bool,
}

impl IntersectionComputation {
    pub fn new(r: &Ray4D, i: &Intersection) -> IntersectionComputation {
        let t = i.t;
        let obj = Rc::clone(&i.what);
        let point = r.position(t);
        let eyev = -r.direction;
        let mut normalv = obj.borrow().normal(point);

        let inside = if normalv.dot(&eyev) < 0.0 {
            normalv = -normalv;
            true
        } else {
            false
        };

        IntersectionComputation { t, obj, point, eyev, normalv, inside }
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
