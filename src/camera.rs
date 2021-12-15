use crate::ray::Ray4D;
use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::world::World;
use crate::canvas::Canvas;
use crate::REFLECTION_RECURSION_DEPTH;

/// A camera record for generating a canvas.
///
/// This record gives a "frame" of the world. Based on camera parameters,
/// different perspectives can be produced.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Camera {
    /// The horizontal size of the resultant canvas.
    pub hsize: usize,

    /// The vertical size of the resultant canvas.
    pub vsize: usize,

    pub half_width: f64,
    pub half_height: f64,
    pub pixel_size: f64,
 
    /// The angle describing "how much" the camera can see.
    pub field_of_view: f64,

    /// A matrix describing how the world should be oriented relative to the
    /// camera (typically a view transformation).
    pub transform: Matrix4D,
}

impl Camera {
    pub fn new(hsize: usize, vsize: usize, field_of_view: f64,
        transform: Matrix4D) -> Camera {
        let half_view = (field_of_view / 2.0).tan();
        let aspect = (hsize as f64) / (vsize as f64);

        let half_width: f64;
        let half_height: f64;

        if aspect >= 1.0 {
            half_width = half_view;
            half_height = half_view / aspect;
        } else {
            half_width = half_view * aspect;
            half_height = half_view;
        }

        let pixel_size = half_width * 2.0 / (hsize as f64);
        Camera {
            hsize,
            vsize,
            half_width,
            half_height,
            pixel_size,
            field_of_view,
            transform,
        }
    }

    pub fn ray_for_pixel(&self, px: usize, py: usize) -> Ray4D {
        // Offsets from the edge of the canvas to the pixel's center
        let xoffset = (px as f64 + 0.5) * self.pixel_size;
        let yoffset = (py as f64 + 0.5) * self.pixel_size;

        // The untransformed coordinates of the pixel in world space
        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        // Using the camera matrix, transform the canvas point and origin,
        // computin the ray's direction vector
        let tr_inv = self.transform.inverse()
            .expect("Camera matrix should have an inverse.");
        
        let pixel = tr_inv * Tuple4D::point(world_x, world_y, -1.0);
        let origin = tr_inv * Tuple4D::point(0.0, 0.0, 0.0);
        let direction = (pixel - origin).normalize();

        Ray4D::new(origin, direction)
    }

    pub fn render(&self, w: &mut World) -> Canvas {
        let mut image = Canvas::new(self.hsize, self.vsize);

        for y in 0..self.vsize {
            for x in 0..self.hsize {
                let ray = self.ray_for_pixel(x, y);
                let color = w.color_at(ray, REFLECTION_RECURSION_DEPTH);
                image.write_pixel(x, y, &color);
            }
        }

        image
    }
}

#[test]
fn ray_through_center() {
    let c = Camera::new(201, 101, std::f64::consts::PI / 2.0,
        Matrix4D::identity());
    let r = c.ray_for_pixel(100, 50);

    assert_eq!(r.origin, Tuple4D::point(0.0, 0.0, 0.0));
    assert_eq!(r.direction, Tuple4D::vector(0.0, 0.0, -1.0));
}

#[test]
fn ray_through_corner() {
    let c = Camera::new(201, 101, std::f64::consts::PI / 2.0,
        Matrix4D::identity());
    let r = c.ray_for_pixel(0, 0);

    assert_eq!(r.origin, Tuple4D::point(0.0, 0.0, 0.0));
    assert_eq!(r.direction, Tuple4D::vector(0.66519, 0.33259, -0.66851));
}

#[test]
fn ray_when_camera_transformed() {
    let c = Camera::new(201, 101, std::f64::consts::PI / 2.0,
        Matrix4D::rotation_y(std::f64::consts::PI / 4.0)
            * Matrix4D::translation(0.0, -2.0, 5.0));
    let r = c.ray_for_pixel(100, 50);

    assert_eq!(r.origin, Tuple4D::point(0.0, 2.0, -5.0));
    assert_eq!(r.direction,
        Tuple4D::vector(2.0f64.sqrt() / 2.0, 0.0, -(2.0f64.sqrt() / 2.0)));
}

#[test]
fn render_world_with_camera() {
    use crate::color::Color;

    let mut w: World = Default::default();
    let mut c = Camera::new(11, 11, std::f64::consts::PI / 2.0,
        Matrix4D::identity());

    let from = Tuple4D::point(0.0, 0.0, -5.0);
    let to = Tuple4D::point(0.0, 0.0, 0.0);
    let up = Tuple4D::vector(0.0, 1.0, 0.0);

    c.transform = Matrix4D::view_transform(from, to, up);
    
    let image = c.render(&mut w);
    assert_eq!(image.read_pixel(5, 5).unwrap(),
        Color::rgb(0.38066, 0.47583, 0.2855));
}
