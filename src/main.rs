use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;
use ray_tracer_challenge::extra::hexagon;

const CANVAS_WIDTH: usize = 960;
const CANVAS_HEIGHT: usize = 540; 

fn main() {
    let mut world = World::empty();

    world.light_source = PointLight::new(
        Color::rgb(0.6, 0.8, 1.0),
        Tuple4D::point(-10.0, 10.0, -10.0),
    );

    let hex = hexagon();
    hex.borrow_mut().transform =
        Matrix4D::rotation_x(2.0 * std::f64::consts::PI / 3.0);

    world.objects = vec![
        hex
    ];

    let mut camera = Camera::new(CANVAS_WIDTH, CANVAS_HEIGHT,
        std::f64::consts::PI / 3.0, Matrix4D::identity());
    camera.transform = Matrix4D::view_transform(
        Tuple4D::point(0.0, 1.5, -5.0),
        Tuple4D::point(0.0, 1.0,  0.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
    ) * Matrix4D::translation(0.0, 0.0, 10.0);

    let canvas = camera.render(&mut world);
    canvas.save("out.ppm").unwrap();
}
