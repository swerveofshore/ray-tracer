use crate::color::Color;
use crate::pattern::Pattern;
use crate::tuple::Tuple4D;
use crate::shape::Shape;

/// A point light.
///
/// A very simple light source. Provides a color and a position where light is
/// produced from.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct PointLight {
    pub intensity: Color,
    pub position: Tuple4D,
}

impl PointLight {
    /// Creates a point light.
    ///
    /// If `position` isn't a point, it is converted to a point automatically.
    pub fn new(intensity: Color, mut position: Tuple4D) -> PointLight {
        if !position.is_point() {
            position.w = 1.0;
        }

        PointLight { intensity, position }
    }
}

/// A material record.
///
/// Materials use attributes from the Phong reflection model; ambient, diffuse,
/// specular and shininess.
#[derive(Clone, Debug, PartialEq)]
pub struct Material {
    /// The inherent color of the material. Ignored if a pattern is specified.
    pub color: Color,

    /// The pattern of a material. See `Pattern` for more information.
    pub pattern: Option<Pattern>,

    /// The ambient light reflected by the material.
    pub ambient: f64,

    /// The diffuse light reflected by the material.
    pub diffuse: f64,

    /// The specular light reflected by the material.
    pub specular: f64,

    /// The shininess of the material (how much white light is shown).
    pub shininess: f64,

    /// The reflectivity of the material (how much light is reflected onto other
    /// objects).
    pub reflective: f64,

    /// The refractive index of the material (degree to which light changes
    /// direction as it passes from one medium to this material).
    pub refractive_index: f64,

    /// The transparency of the material.
    pub transparency: f64,
}

/// Defaults for a material.
///
/// By default, a material is white with no pattern. It has some ambient light,
/// reflects almost all diffuse and specular light, and is relatively shiny.
///
/// The default material is non-reflective, has the refractive index of a
/// vacuum (refraction is largely irrelevant). The material is opaque.
///
/// # Examples
///
/// Creating a default material:
///
/// ```
/// # use ray_tracer::color::Color;
/// # use ray_tracer::light::Material;
/// let default_material: Material = Default::default();
/// let actual_defaults = Material {
///     color: Color::white(),
///     pattern: None,
///
///     ambient: 0.1,
///     diffuse: 0.9,
///     specular: 0.9,
///     shininess: 200.0,
///
///     reflective: 0.0,
///     refractive_index: 1.0,
///     transparency: 0.0
/// };
/// assert_eq!(default_material, actual_defaults);
/// ```
impl Default for Material {
    fn default() -> Material {
        Material {
            color: Color::rgb(1.0, 1.0, 1.0),
            pattern: None,

            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,

            reflective: 0.0,
            refractive_index: 1.0,
            transparency: 0.0,
        }
    }
}

/// Calculate the lighting of a pixel in an environment.
///
/// Inline comments are provided in the source for this function. I believe they
/// best explain how this function operates.
///
/// Effectively, this function takes a material, a single light, a point, the
/// eye vector and the normal vector, and calculates how the light looks from
/// the eye. Position is irrelevant, bar the angle around `point`.
///
/// If this point is in a shadow (parameter `in_shadow`), only ambient light is
/// used.
pub fn lighting(m: &Material, obj: &Shape, light: &PointLight,
    point: Tuple4D, eyev: Tuple4D, normalv: Tuple4D, in_shadow: bool) -> Color {
    // If Material m has some pattern, use that for color.
    let color = if let Some(ref pat) = m.pattern {
        pat.pattern_at_object(obj, point) 
    } else {
        obj.material().color
    };

    // Combine surface color with light's color.
    let effective_color = color * light.intensity;

    // Find direction to light source.
    let lightv = (light.position - point).normalize();

    // Compute ambient light.
    let ambient = effective_color * m.ambient;

    // If the point is in a shadow, only calculate ambient light.
    if in_shadow {
        return ambient;
    }

    // Declare diffuse and specular variables for calculating light.
    let diffuse;
    let specular;

    // For the side of the surface with no light, use only ambient light.
    let light_dot_normal = lightv.dot(&normalv);
    if light_dot_normal < 0.0 {
        diffuse = Color::black();
        specular = Color::black();
    } else {
        diffuse = effective_color * m.diffuse * light_dot_normal;

        let reflectv = (-lightv).reflect(&normalv);
        let reflect_dot_eye = reflectv.dot(&eyev);

        // If no specular reflection, set the specular light to black.
        if reflect_dot_eye <= 0.0 {
            specular = Color::black();
        } else {
            // Otherwise, calculate the shininess factor and apply.
            let factor = reflect_dot_eye.powf(m.shininess);
            specular = light.intensity * m.specular * factor;
        }
    }

    ambient + diffuse + specular
}

#[test]
fn lighting_with_eye_between_light_and_surface() {
    use crate::shape::Shape;

    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);
    let mut s = Shape::sphere();
    s.material = m;

    let eyev = Tuple4D::vector(0.0, 0.0, -1.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, -10.0),
    );

    let res = lighting(
        &s.material, &s, &light, position, eyev, normalv, false
    );
    assert_eq!(res, Color::rgb(1.9, 1.9, 1.9));
}

#[test]
fn lighting_with_eye_between_light_and_surface_offset_45() {
    use crate::shape::Shape;

    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);
    let mut s = Shape::sphere();
    s.material = m;

    let eyev = Tuple4D::vector(0.0, 2.0f64.sqrt() / 2.0, 2.0f64.sqrt() / 2.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, -10.0),
    );

    let res = lighting(
        &s.material, &s, &light, position, eyev, normalv, false
    );
    assert_eq!(res, Color::rgb(1.0, 1.0, 1.0));
}

#[test]
fn lighting_with_eye_opposite_from_surface_offset_45() {
    use crate::shape::Shape;

    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);
    let mut s = Shape::sphere();
    s.material = m;

    let eyev = Tuple4D::vector(0.0, 0.0, -1.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 10.0, -10.0),
    );

    let res = lighting(
        &s.material, &s, &light, position, eyev, normalv, false
    );
    assert_eq!(res, Color::rgb(0.7364, 0.7364, 0.7364));
}

#[test]
fn lighting_with_eye_opposite_from_surface_in_reflection() {
    use crate::shape::Shape;

    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);
    let mut s = Shape::sphere();
    s.material = m;

    let eyev = Tuple4D::vector(0., -(2.0f64.sqrt())/2., -(2.0f64.sqrt())/2.);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 10.0, -10.0),
    );

    let res = lighting(
        &s.material, &s, &light, position, eyev, normalv, false
    );
    assert_eq!(res, Color::rgb(1.6364, 1.6364, 1.6364));
}

#[test]
fn lighting_with_eye_across_surface_from_light() {
    use crate::shape::Shape;

    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);
    let mut s = Shape::sphere();
    s.material = m;

    let eyev = Tuple4D::vector(0.0, 0.0, -1.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, 10.0),
    );

    let res = lighting(
        &s.material, &s, &light, position, eyev, normalv, false
    );
    assert_eq!(res, Color::rgb(0.1, 0.1, 0.1));   
}

#[test]
fn lighting_with_stripe_pattern() {
    use crate::shape::Shape;

    let m = Material {
        color: Color::rgb(0.5, 0.5, 0.5),
        pattern: Some(
            Pattern::stripe(Color::white(),Color::black())
        ),

        // Note that ONLY ambient light is included, as the color of ambient
        // light is mostly predictable
        ambient: 1.0,
        diffuse: 0.0,
        specular: 0.0,
        shininess: 0.0,
        reflective: 0.0,

        ..Default::default()
    };

    let mut s = Shape::sphere();
    s.material = m;

    let eyev = Tuple4D::vector(0.0, 0.0, -1.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::white(), Tuple4D::point(0.0, 0.0, -10.0)
    );

    assert_eq!(
        Color::white(),
        lighting(&s.material, &s, &light, Tuple4D::point(0.9, 0.0, 0.0),
            eyev, normalv, false)
    );

    assert_eq!(
        Color::black(),
        lighting(&s.material, &s, &light, Tuple4D::point(1.1, 0.0, 0.0),
            eyev, normalv, false)
    );
}
