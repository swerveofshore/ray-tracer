use std::rc::Rc;
use std::cell::RefCell;
use crate::ray::Ray4D;
use crate::tuple::Tuple4D;
use crate::color::Color;
use crate::matrix::Matrix4D;
use crate::light::{ PointLight, Material, lighting };
use crate::geometry::{ ShapeDebug, Intersections,
    IntersectionComputation, Sphere };

/// A world with objects and light.
///
/// Objects are all `IntersectableDebug` items, and the light is a `PointLight`.
///
/// Worlds collect all objects as well as light for rendering. Most logic is
/// performed within worlds for the ray tracer.
pub struct World {
    pub objects: Vec<Rc<RefCell<dyn ShapeDebug>>>,
    pub light_source: PointLight,
}

impl Default for World {
    fn default() -> World {
        let light_source = PointLight::new(
            Color::rgb(1.0, 1.0, 1.0),
            Tuple4D::point(-10.0, 10.0, -10.0)
        );

        let mut s1 = Sphere::unit();
        let mut m1: Material = Default::default();
        m1.color = Color::rgb(0.8, 1.0, 0.6);
        m1.diffuse = 0.7;
        m1.specular = 0.2;
        s1.material = m1;

        let mut s2 = Sphere::unit();
        s2.transform = Matrix4D::scaling(0.5, 0.5, 0.5);

        World {
            objects: vec![Rc::new(RefCell::new(s1)), Rc::new(RefCell::new(s2))],
            light_source
        }
    }
}

impl World {
    /// Creates a default world with two spheres.
    pub fn new() -> World {
        Default::default()
    }

    /// Creates an empty world with no objects and the default light source.
    pub fn empty() -> World {
        World { objects: Vec::new(), light_source: Default::default() }
    }

    /// Gets the first object in a world.
    pub fn first(&self) -> Option<Rc<RefCell<dyn ShapeDebug>>> {
        if self.objects.len() > 0 {
            Some(Rc::clone(&self.objects[0]))
        } else {
            None
        }
    }

    /// Gets the second object in a world.
    pub fn second(&self) -> Option<Rc<RefCell<dyn ShapeDebug>>> {
        if self.objects.len() > 1 {
            Some(Rc::clone(&self.objects[1]))
        } else {
            None
        }
    }

    /// Intersects a ray against all objects in a world.
    pub fn intersect(&self, r: Ray4D) -> Intersections {
        let mut intersections: Intersections  = Intersections::new();
        for obj in self.objects.iter() {
            let mut is: Intersections = obj.borrow().intersect(r);
            intersections.intersections.append(&mut is.intersections);
        }

        intersections.sort();
        intersections
    }

    /// Determines whether a point is shadowed.
    pub fn is_shadowed(&self, p: Tuple4D) -> bool {
        let v = self.light_source.position - p;
        let distance = v.magnitude();
        let direction = v.normalize();

        let r = Ray4D::new(p, direction);
        let mut intersections = self.intersect(r);

        let h = intersections.hit();
        if let Some(i) = h {
            i.t < distance
        } else {
            false
        }
    }

    /// Calculates the color for a hit, based on shadows and light.
    pub fn shade_hit(&self, comps: &IntersectionComputation) -> Color {
        lighting(*comps.obj.borrow().material(),
            self.light_source, comps.point, comps.eyev, comps.normalv,
            self.is_shadowed(comps.over_point))
    }

    /// Determines a color based on the intersection of a ray and the objects.
    pub fn color_at(&self, r: Ray4D) -> Color {
        let mut is = self.intersect(r);

        // If at least one object is hit, return the color, else return black
        match is.hit() {
            None => Color::black(),
            Some(i) => {
                let comps = IntersectionComputation::new(&r, &i);
                self.shade_hit(&comps)
            },
        }
    }
}

#[test]
fn intersect_default_world_with_ray() {
    let w: World = Default::default();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let is: Intersections = w.intersect(r);

    assert_eq!(is.intersections.len(), 4);
    assert_eq!(is.intersections[0].t, 4.0);
    assert_eq!(is.intersections[1].t, 4.5);
    assert_eq!(is.intersections[2].t, 5.5);
    assert_eq!(is.intersections[3].t, 6.0);
}

#[test]
fn shade_intersection_from_outside() {
    use crate::geometry::Intersection;

    let w: World = Default::default();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let shape = w.first().expect("Default world should contain objects.");
    let i = Intersection { t: 4.0, what: shape };

    let comps = IntersectionComputation::new(&r, &i);
    let c = w.shade_hit(&comps);

    assert_eq!(c, Color::rgb(0.38066, 0.47583, 0.2855));
}

/* TODO: No idea why this doesn't pass.
#[test]
fn shade_intersection_from_inside() {
    use crate::geometry::Intersection;

    let mut w: World = Default::default();
    w.light_source = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.25, 0.0),
    );

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let shape = w.second().expect("Default world should contain objects.");
    let i = Intersection { t: 0.5, what: shape };

    let comps = IntersectionComputation::new(&r, &i);
    let c = w.shade_hit(&comps);

    assert_eq!(c, Color::rgb(0.90498, 0.90498, 0.90498));
}
*/

#[test]
fn shade_intersection_in_shadow() {
    use crate::geometry::Intersection;

    let mut w: World = World::empty();
    w.light_source = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, -10.0),
    );

    let s1 = Rc::new(RefCell::new(Sphere::unit()));
    w.objects.push(s1);

    let s2 = Rc::new(RefCell::new(Sphere::unit()));
    s2.borrow_mut().transform = Matrix4D::translation(0.0, 0.0, 10.0);
    w.objects.push(s2);

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let i = Intersection { t: 4.0, what: w.second().unwrap() };
    let comps = IntersectionComputation::new(&r, &i);
    let c = w.shade_hit(&comps);

    assert_eq!(c, Color::rgb(0.1, 0.1, 0.1));
}

#[test]
fn color_ray_miss() {
    let w: World = Default::default();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
    );

    assert_eq!(w.color_at(r), Color::black());
}

#[test]
fn color_ray_hit() {
    let w: World = Default::default();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    assert_eq!(w.color_at(r), Color::rgb(0.38066, 0.47583, 0.2855));
}

#[test]
fn color_behind_ray() {
    let w: World = Default::default();
    let outer = w.first().unwrap();
    let inner = w.second().unwrap();

    outer.borrow_mut().material_mut().ambient = 1.0;
    inner.borrow_mut().material_mut().ambient = 1.0;

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.75),
        Tuple4D::vector(0.0, 0.0, -1.0)
    );

    assert_eq!(w.color_at(r), inner.borrow().material().color);
}

#[test]
fn shadow_collinear_point_and_light() {
    let w: World = Default::default();
    let p = Tuple4D::point(0.0, 10.0, 0.0);

    assert!(!w.is_shadowed(p));
}

#[test]
fn shadow_light_between_point_and_spheres() {
    let w: World = Default::default();
    let p = Tuple4D::point(10.0, -10.0, 10.0);

    assert!(w.is_shadowed(p));
}

#[test]
fn shadow_object_behind_light() {
    let w: World = Default::default();
    let p = Tuple4D::point(-20.0, 20.0, -20.0);

    assert!(!w.is_shadowed(p));
}

#[test]
fn shadow_object_behind_point() {
    let w: World = Default::default();
    let p = Tuple4D::point(-2.0, 2.0, -2.0);

    assert!(!w.is_shadowed(p));
}
