use crate::feq;
use crate::ray::Ray4D;
use crate::tuple::Tuple4D;
use crate::color::Color;
use crate::matrix::Matrix4D;
use crate::light::{ PointLight, Material, lighting };
use crate::geometry::{ ShapeDebug, Sphere };
use crate::intersect::{ Intersections, IntersectionComputation };
use crate::shape::intersect;

/// A world with objects and light.
///
/// Objects are all `IntersectableDebug` items, and the light is a `PointLight`.
///
/// Worlds collect all objects as well as light for rendering. Most logic is
/// performed within worlds for the ray tracer.
pub struct World {
    pub objects: Vec<Box<dyn ShapeDebug>>,
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
            objects: vec![Box::new(s1), Box::new(s2)],
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

    /// Gets a reference to the first object in a world.
    pub fn first(&self) -> Option<& dyn ShapeDebug> {
        if self.objects.len() > 0 {
            Some(&*self.objects[0])
        } else {
            None
        }
    }

    /// Gets a mutable reference to the first object in a world.
    pub fn first_mut(&mut self) -> Option<&mut dyn ShapeDebug> {
        if self.objects.len() > 0 {
            Some(&mut *self.objects[0])
        } else {
            None
        }
    }

    /// Gets a reference to the second object in a world.
    pub fn second(&self) -> Option<& dyn ShapeDebug> {
        if self.objects.len() > 1 {
            Some(&*self.objects[1])
        } else {
            None
        }
    }

    /// Gets a mutable reference to the second objet in a world.
    pub fn second_mut(&mut self) -> Option<&mut dyn ShapeDebug> {
        if self.objects.len() > 1 {
            Some(&mut *self.objects[1])
        } else {
            None
        }
    }

    /// Intersects a ray against all objects in a world.
    ///
    /// The `World` needs to be mutable because it owns the objects it contains.
    /// Function `intersect` has the opportunity to change the `saved_ray`
    /// associated with each object.
    pub fn intersect(&self, r: Ray4D) -> Intersections {
        let mut intersections: Intersections = Intersections::new();
        for obj in self.objects.iter() {
            let mut is: Intersections = intersect(&**obj, r);
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

    pub fn reflected_color(&self, comps: &IntersectionComputation,
        remaining: usize) -> Color {

        // If there are no more recursive calls remaining, exit.
        if remaining == 0 {
            return Color::black();
        }

        // If the object is non-reflective, return black (no reflected color).
        if feq(comps.obj.material().reflective, 0.0) {
            return Color::black();
        }

        let reflect_ray = Ray4D::new(comps.over_point, comps.reflectv);
        let color = self.color_at(reflect_ray, remaining - 1);

        color * comps.obj.material().reflective
    }

    /// Returns the refracted color of an object, given sufficient transparency.
    ///
    /// An object is considered "opaque" if its material's transparency property
    /// is zero. In that case, no refraction occurs and black is returned.
    ///
    /// Sometimes refraction doesn't "escape" the object, particularly under
    /// "total internal reflection." When light enters another medium at a
    /// sufficiently acute angle, and the new medium has a lower refractive
    /// index than the old, the ray may reflect beneath the surface eternally.
    ///
    /// This is due to Snell's law. For reference, the law states:
    ///
    /// ```latex
    /// \frac{\sin \theta_i}{\sin \theta_t} = \frac{\nu_2}{\nu_1}
    /// ```
    ///
    /// In other words, the angle of the ray as it exits the surface, divided
    /// the angle of the ray as it leaves/re-enters the surface, is equal to the
    /// outer refraction constant divided by the first. This is used to
    /// determine whether "total internal reflection" is occurring.
    pub fn refracted_color(&self, comps: &IntersectionComputation,
        remaining: usize) -> Color {

        // If there are no more recursive calls remaining, exit.
        if remaining == 0 {
            return Color::black();
        }

        // If the object is opaque, return black (no refracted color).
        if feq(comps.obj.material().transparency, 0.0) {
            return Color::black()
        }

        let n_ratio = comps.n1 / comps.n2;
        let cos_i = comps.eyev.dot(&comps.normalv);
        let sin2_t = n_ratio.powi(2) * (1.0 - cos_i.powi(2));

        // If total internal reflection is occurring, return black.
        if sin2_t > 1.0 {
            return Color::black();
        }

        let cos_t = (1.0 - sin2_t).sqrt();

        // Direction of the refracted ray.
        let direction = comps.normalv * (n_ratio * cos_i - cos_t)
            - comps.eyev * n_ratio;

        // The refracted ray.
        let refracted_ray = Ray4D::new(comps.under_point, direction);

        // Color of the refracted ray, multiplied by transparency for opacity.
        self.color_at(refracted_ray, remaining - 1)
            * comps.obj.material().transparency
    }

    /// Calculates the color for a hit, based on shadows and light.
    pub fn shade_hit(&self, comps: &IntersectionComputation, remaining: usize)
        -> Color {
        // Color calculated from light on a surface
        let surface = lighting(
            *comps.obj.material(), comps.obj, self.light_source,
             comps.over_point, comps.eyev, comps.normalv,
             self.is_shadowed(comps.over_point)
        );

        // Color calculated from light reflecting off of a surface
        let reflected = self.reflected_color(&comps, remaining);

        // Color calculated from light going through a surface
        let refracted = self.refracted_color(&comps, remaining);

        // If a material is reflective and transparent, use Schlick approx.
        let m = comps.obj.material();
        if m.reflective > 0.0 && m.transparency > 0.0 {
            let reflectance = comps.schlick();

            // Use the approximation.
            surface + reflected * reflectance + refracted * (1.0 - reflectance)
        } else {
            // Add together all calculated colors.
            surface + reflected + refracted
        }
    }

    /// Determines a color based on the intersection of a ray and the objects.
    pub fn color_at(&self, r: Ray4D, remaining: usize) -> Color {
        let mut is = self.intersect(r);
        let hit = is.hit();

        // If at least one object is hit, return the color, else return black
        match hit {
            None => Color::black(),
            Some(i) => {
                let comps = IntersectionComputation::new(&r, &i, Some(&is));
                self.shade_hit(&comps, remaining)
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
    use crate::consts::REFLECTION_RECURSION_DEPTH;
    use crate::intersect::Intersection;

    let w: World = Default::default();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let shape = w.first().expect("Default world should contain objects.");
    let i = Intersection { t: 4.0, what: shape };

    let comps = IntersectionComputation::new(&r, &i, None);
    let c = w.shade_hit(&comps, REFLECTION_RECURSION_DEPTH);

    assert_eq!(c, Color::rgb(0.38066, 0.47583, 0.2855));
}

#[test]
fn shade_intersection_from_inside() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;
    use crate::intersect::Intersection;

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

    let comps = IntersectionComputation::new(&r, &i, None);
    let c = w.shade_hit(&comps, REFLECTION_RECURSION_DEPTH);

    assert_eq!(c, Color::rgb(0.90498, 0.90498, 0.90498));
}

#[test]
fn shade_intersection_in_shadow() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;
    use crate::intersect::Intersection;

    let mut w: World = World::empty();
    w.light_source = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, -10.0),
    );

    let s1 = Sphere::unit();
    w.objects.push(Box::new(s1));

    let mut s2 = Sphere::unit();
    s2.transform = Matrix4D::translation(0.0, 0.0, 10.0);
    w.objects.push(Box::new(s2));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    let i = Intersection { t: 4.0, what: w.second().unwrap() };
    let comps = IntersectionComputation::new(&r, &i, None);
    let c = w.shade_hit(&comps, REFLECTION_RECURSION_DEPTH);

    assert_eq!(c, Color::rgb(0.1, 0.1, 0.1));
}

#[test]
fn color_ray_miss() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;

    let w: World = Default::default();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
    );

    assert_eq!(w.color_at(r, REFLECTION_RECURSION_DEPTH), Color::black());
}

#[test]
fn color_ray_hit() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;

    let w: World = Default::default();
    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0),
    );

    assert_eq!(
        w.color_at(r, REFLECTION_RECURSION_DEPTH),
        Color::rgb(0.38066, 0.47583, 0.2855)
    );
}

#[test]
fn color_behind_ray() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;

    let mut w: World = Default::default();
    {
        let outer = w.first_mut().unwrap();
        outer.material_mut().ambient = 1.0;
    }
    {
        let inner_mut = w.second_mut().unwrap();
        inner_mut.material_mut().ambient = 1.0;
    }

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.75),
        Tuple4D::vector(0.0, 0.0, -1.0)
    );

    let inner_color = {
        let inner = w.second().unwrap();
        inner.material().color
    };
    assert_eq!(w.color_at(r, REFLECTION_RECURSION_DEPTH), inner_color);
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

#[test]
fn reflected_color_for_nonreflective_material() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;
    use crate::intersect::Intersection;

    let mut w: World = Default::default();

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    {
        let s = w.second_mut().unwrap();
        s.material_mut().ambient = 1.0;
    }

    let s = w.second().unwrap();
    let i = Intersection { t: 1.0, what: s };
    let comps = IntersectionComputation::new(&r, &i, None);

    assert_eq!(
        Color::black(),
        w.reflected_color(&comps, REFLECTION_RECURSION_DEPTH)
    );
}

#[test]
fn reflected_color_for_reflective_material() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;
    use crate::intersect::Intersection;
    use crate::geometry::Plane;

    let mut s = Plane::new();
    s.material.reflective = 0.5;
    s.transform = Matrix4D::translation(0.0, -1.0, 0.0);

    let mut w: World = Default::default();
    w.objects.push(Box::new(s));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -3.0),
        Tuple4D::vector(0.0, -(2.0f64.sqrt()) / 2.0, (2.0f64.sqrt()) / 2.0)
    );

    let i = Intersection { t: 2.0f64.sqrt(), what: &*w.objects[2] };
    let comps = IntersectionComputation::new(&r, &i, None);

    assert_eq!(
        Color::rgb(0.19032, 0.2379, 0.14274),
        w.reflected_color(&comps, REFLECTION_RECURSION_DEPTH)
    );
}

#[test]
fn shade_hit_with_reflective_material() {
    use crate::consts::REFLECTION_RECURSION_DEPTH;
    use crate::intersect::Intersection;
    use crate::geometry::Plane;

    let mut s = Plane::new();
    s.material.reflective = 0.5;
    s.transform = Matrix4D::translation(0.0, -1.0, 0.0);

    let mut w: World = Default::default();
    w.objects.push(Box::new(s));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -3.0),
        Tuple4D::vector(0.0, -(2.0f64.sqrt()) / 2.0, (2.0f64.sqrt()) / 2.0)
    );

    let i = Intersection { t: 2.0f64.sqrt(), what: &*w.objects[2] };
    let comps = IntersectionComputation::new(&r, &i, None);

    assert_eq!(
        Color::rgb(0.87677, 0.92436, 0.82918),
        w.shade_hit(&comps, REFLECTION_RECURSION_DEPTH)
    );
}

#[test]
fn refracted_color_on_opaque_material() {
    use crate::consts::REFRACTION_RECURSION_DEPTH;
    use crate::intersect::{ Intersection, Intersections };

    let w: World = Default::default();
    let s = w.first().unwrap();

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: 4.0, what: s },
            Intersection { t: 6.0, what: s },
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[0], Some(&is)
    );

    assert_eq!(
        Color::black(),
        w.refracted_color(&comps, REFRACTION_RECURSION_DEPTH)
    );
}

#[test]
fn refracted_color_at_max_recursion_depth() {
    use crate::intersect::{ Intersection, Intersections };

    let mut w: World = Default::default();
    {
        let s = w.first_mut().unwrap();
        s.material_mut().transparency = 1.0;
        s.material_mut().refractive_index = 1.5;
    }
    let s = w.first().unwrap();

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -5.0),
        Tuple4D::vector(0.0, 0.0, 1.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: 4.0, what: s },
            Intersection { t: 6.0, what: s },
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[0], Some(&is)
    );

    assert_eq!(
        Color::black(),
        w.refracted_color(&comps, 0)
    );
}

#[test]
fn refracted_color_under_total_internal_reflection() {
    use crate::consts::REFRACTION_RECURSION_DEPTH;
    use crate::intersect::Intersection;

    let mut w: World = Default::default();
    {
        let s = w.first_mut().unwrap();
        s.material_mut().transparency = 1.0;
        s.material_mut().refractive_index = 1.5;
    }
    let s = w.first().unwrap();

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 2.0f64.sqrt() / 2.0),
        Tuple4D::vector(0.0, 1.0, 0.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: -(2.0f64.sqrt()) / 2.0, what: s },
            Intersection { t: 2.0f64.sqrt() / 2.0, what: s },
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[1], Some(&is)
    );

    assert_eq!(
        Color::black(),
        w.refracted_color(&comps, REFRACTION_RECURSION_DEPTH)
    );
}

#[test]
fn refracted_color_with_refracted_ray() {
    use crate::intersect::Intersection;
    use crate::pattern::Pattern;

    let mut w: World = Default::default();
    {
        let a = w.first_mut().unwrap();
        a.material_mut().ambient = 1.0;
        a.material_mut().pattern = Some(Pattern::point());
    }
    {
        let b = w.second_mut().unwrap();
        b.material_mut().transparency = 1.0;
        b.material_mut().refractive_index = 1.5;
    }

    let a = w.first().unwrap();
    let b = w.second().unwrap();

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, 0.1),
        Tuple4D::vector(0.0, 1.0, 0.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: -0.9899, what: a },
            Intersection { t: -0.4899, what: b },
            Intersection { t:  0.4899, what: b },
            Intersection { t:  0.9899, what: a },
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[2], Some(&is)
    );

    assert_eq!(
        Color::rgb(0.0, 0.99888, 0.04725),
        w.refracted_color(&comps, 5)
    );
}

#[test]
fn shade_hit_with_transparent_material() {
    use crate::intersect::Intersection;
    use crate::geometry::Plane;

    let mut w: World = Default::default();

    let mut floor = Plane::new();
    floor.transform = Matrix4D::translation(0.0, -1.0, 0.0);
    floor.material.transparency = 0.5;
    floor.material.refractive_index = 1.5;

    let mut ball = Sphere::unit();
    ball.material.color = Color::red();
    ball.material.ambient = 0.5;
    ball.transform = Matrix4D::translation(0.0, -3.5, -0.5);

    w.objects.push(Box::new(floor));
    w.objects.push(Box::new(ball));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -3.0),
        Tuple4D::vector(0.0, -(2.0f64.sqrt()) / 2.0, 2.0f64.sqrt() / 2.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: 2.0f64.sqrt(), what: &*w.objects[2] }
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[0], Some(&is)
    );

    assert_eq!(
        Color::rgb(0.93642, 0.68642, 0.68642),
        w.shade_hit(&comps, 5)
    );
}

#[test]
fn shade_hit_with_reflective_transparent_material() {
    use crate::intersect::Intersection;
    use crate::geometry::Plane;

    let mut w: World = Default::default();

    let mut floor = Plane::new();
    floor.transform = Matrix4D::translation(0.0, -1.0, 0.0);
    floor.material.reflective = 0.5;
    floor.material.transparency = 0.5;
    floor.material.refractive_index = 1.5;

    let mut ball = Sphere::unit();
    ball.material.color = Color::red();
    ball.material.ambient = 0.5;
    ball.transform = Matrix4D::translation(0.0, -3.5, -0.5);

    w.objects.push(Box::new(floor));
    w.objects.push(Box::new(ball));

    let r = Ray4D::new(
        Tuple4D::point(0.0, 0.0, -3.0),
        Tuple4D::vector(0.0, -(2.0f64.sqrt()) / 2.0, 2.0f64.sqrt() / 2.0)
    );

    let is = Intersections {
        intersections: vec![
            Intersection { t: 2.0f64.sqrt(), what: &*w.objects[2] }
        ]
    };

    let comps = IntersectionComputation::new(
        &r, &is.intersections[0], Some(&is)
    );

    assert_eq!(
        Color::rgb(0.93391, 0.69643, 0.69243),
        w.shade_hit(&comps, 5)
    );
}
