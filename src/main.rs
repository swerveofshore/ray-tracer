use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::geometry::{ Sphere, Plane };
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;

const CANVAS_WIDTH: usize = 3840;
const CANVAS_HEIGHT: usize = 2160;

fn main() {
    let mut floor = Plane::new();
    floor.transform = Matrix4D::scaling(10.0, 0.01, 10.0);
    floor.material = Default::default();
    floor.material.color = Color::rgb(0.5, 0.5, 0.5);
    floor.material.specular = 0.0;

    let mut middle = Sphere::unit();
    middle.transform = Matrix4D::translation(-0.5, 1.0, 2.0);
    middle.material = Default::default();
    middle.material.color = Color::rgb(1.0, 0.4666, 0.2666);
    middle.material.diffuse = 0.7;
    middle.material.specular = 0.3;

    let mut right = Sphere::unit();
    right.transform = Matrix4D::translation(1.5, 0.5, -0.5)
        * Matrix4D::scaling(0.5, 0.5, 0.5);
    right.material = Default::default();
    right.material.color = Color::rgb(1.0, 0.6666, 0.2666);
    right.material.diffuse = 0.7;
    right.material.specular = 0.3;

    let mut left = Sphere::unit();
    left.transform = Matrix4D::translation(-1.5, 0.33, -0.75)
        * Matrix4D::scaling(0.33, 0.33, 0.33);
    left.material = Default::default();
    left.material.color = Color::rgb(0.8666, 0.2, 0.2);
    left.material.diffuse = 0.7;
    left.material.specular = 0.3;

    let mut world = World::empty();
    world.light_source = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(-10.0, 10.0, -10.0),
    );
    world.objects = vec![
        Box::new(floor),
        Box::new(middle),
        Box::new(right),
        Box::new(left),
    ];

    let mut camera = Camera::new(CANVAS_WIDTH, CANVAS_HEIGHT,
        std::f64::consts::PI / 3.0, Matrix4D::identity());
    camera.transform = Matrix4D::view_transform(
        Tuple4D::point(0.0, 1.5, -5.0),
        Tuple4D::point(0.0, 1.0,  0.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
    );

    let canvas = camera.render(&mut world);
    canvas.save("out.ppm").unwrap();
}
