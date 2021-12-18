use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::geometry::{ Sphere, Plane, Cube, Cylinder, Cone };
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::pattern::Pattern;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;
const CANVAS_WIDTH: usize = 960 * 4;
const CANVAS_HEIGHT: usize = 540 * 4; 

fn main() {
    let mut floor = Plane::new();
    floor.transform = Matrix4D::scaling(10.0, 0.01, 10.0);
    floor.material = Default::default();
    floor.material.color = Color::rgb(0.5, 0.5, 0.5);
    floor.material.specular = 0.0;
    floor.material.pattern = Some(Pattern::checker(
        Color::white(),
        Color::black()
    ));
    floor.material.pattern.get_or_insert(floor.material.pattern.unwrap())
        .transform = Matrix4D::scaling(0.1, 0.1, 0.1);
    floor.material.reflective = 0.5;

    let mut middle = Sphere::unit();
    middle.transform = Matrix4D::translation(-0.5, 1.0, 2.0);
    middle.material = Default::default();
    middle.material.color = Color::rgb(1.0, 0.4666, 0.2666);
    middle.material.diffuse = 0.7;
    middle.material.specular = 0.3;
    middle.material.transparency = 0.5;
    middle.material.reflective = 0.5;
    middle.material.refractive_index = 1.5;

    /*
    middle.material.pattern = Some(Pattern::checker(
        Color::blue(),
        Color::green()
    ));
    */

    let mut right = Cone::unit();
    right.minimum = 0.0;
    right.maximum = 3.0;
    right.closed = true;
    right.transform = Matrix4D::translation(1.5, 2.5, -0.5)
        * Matrix4D::scaling(0.25, 0.25, 0.25)
        // * Matrix4D::rotation_x(std::f64::consts::PI / 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_y(std::f64::consts::PI / 8.0);
    right.material = Default::default();
    right.material.color = Color::rgb(1.0, 0.6666, 0.2666);
    right.material.diffuse = 0.7;
    right.material.specular = 0.3;
    right.material.reflective = 0.3;

    let mut left = Cube::unit();
    left.transform = Matrix4D::translation(-1.5, 1.0, -0.75)
        * Matrix4D::scaling(0.33, 0.33, 0.33)
        * Matrix4D::rotation_x(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_y(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 4.0);
    left.material = Default::default();
    left.material.color = Color::rgb(0.8666, 0.2, 0.2);
    left.material.diffuse = 0.7;
    left.material.specular = 0.3;
    left.material.transparency = 0.5;
    left.material.refractive_index = 1.1;

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
    ) * Matrix4D::translation(0.0, 0.0, 10.0);

    let canvas = camera.render(&mut world);
    canvas.save("out.ppm").unwrap(); }
