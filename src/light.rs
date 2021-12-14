use crate::color::Color;
use crate::tuple::Tuple4D;

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
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
}

impl Default for Material {
    fn default() -> Material {
        Material {
            color: Color::rgb(1.0, 1.0, 1.0),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0
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
pub fn lighting(m: &Material, light: PointLight, point: Tuple4D,
    eyev: Tuple4D, normalv: Tuple4D, in_shadow: bool) -> Color {
    // Combine surface color with light's color
    let effective_color = m.color * light.intensity;

    // Find direction to light source
    let lightv = (light.position - point).normalize();

    // Compute ambient light
    let ambient = effective_color * m.ambient;

    // If the point is in a shadow, only calculate ambient light
    if in_shadow {
        return ambient;
    }

    // Default to black for locations without light
    let mut diffuse = Color::black();
    let mut specular = Color::black();

    // Find which side of the surface the light is on, act accordingly
    let light_dot_normal = lightv.dot(&normalv);
    if light_dot_normal >= 0.0 {
        diffuse = effective_color * m.diffuse * light_dot_normal;

        let reflectv = (-lightv).reflect(&normalv);
        let reflect_dot_eye = reflectv.dot(&eyev);

        // If no specular reflection, set the specular light to black 
        if reflect_dot_eye <= 0.0 {
            specular = Color::black();
        } else {
            // Otherwise, calculate the shininess factor and apply
            let factor = reflect_dot_eye.powf(m.shininess);
            specular = light.intensity * m.specular * factor;
        }
    }

    ambient + diffuse + specular
}

#[test]
fn eye_between_light_and_surface() {
    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);

    let eyev = Tuple4D::vector(0.0, 0.0, -1.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, -10.0),
    );

    let res = lighting(&m, light, position, eyev, normalv, false);
    assert_eq!(res, Color::rgb(1.9, 1.9, 1.9));
}

#[test]
fn eye_between_light_and_surface_offset_45() {
    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);

    let eyev = Tuple4D::vector(0.0, 2.0f64.sqrt() / 2.0, 2.0f64.sqrt() / 2.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, -10.0),
    );

    let res = lighting(&m, light, position, eyev, normalv, false);
    assert_eq!(res, Color::rgb(1.0, 1.0, 1.0));
}

#[test]
fn eye_opposite_from_surface_offset_45() {
    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);

    let eyev = Tuple4D::vector(0.0, 0.0, -1.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 10.0, -10.0),
    );

    let res = lighting(&m, light, position, eyev, normalv, false);
    assert_eq!(res, Color::rgb(0.7364, 0.7364, 0.7364));
}

#[test]
fn eye_opposite_from_surface_in_reflection() {
    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);

    let eyev = Tuple4D::vector(0., -(2.0f64.sqrt())/2., -(2.0f64.sqrt())/2.);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 10.0, -10.0),
    );

    let res = lighting(&m, light, position, eyev, normalv, false);
    assert_eq!(res, Color::rgb(1.6364, 1.6364, 1.6364));
}

#[test]
fn eye_across_surface_from_light() {
    let m: Material = Default::default();
    let position = Tuple4D::point(0.0, 0.0, 0.0);

    let eyev = Tuple4D::vector(0.0, 0.0, -1.0);
    let normalv = Tuple4D::vector(0.0, 0.0, -1.0);
    let light = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(0.0, 0.0, 10.0),
    );

    let res = lighting(&m, light, position, eyev, normalv, false);
    assert_eq!(res, Color::rgb(0.1, 0.1, 0.1));   
}
