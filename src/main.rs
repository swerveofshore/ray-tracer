use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::shape::Shape;
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::pattern::Pattern;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;

const CANVAS_WIDTH: usize = 960 * 4;
const CANVAS_HEIGHT: usize = 540 * 4; 

fn main() {
    let mut floor = Shape::plane();
    *floor.transform_mut() = Matrix4D::scaling(10.0, 0.01, 10.0);
    *floor.material_mut() = Default::default();
    floor.material_mut().color = Color::rgb(0.5, 0.5, 0.5);
    floor.material_mut().specular = 0.0;
    floor.material_mut().pattern = Some(Pattern::checker(
        Color::white(),
        Color::black()
    ));
    floor.material_mut().pattern.iter_mut().next().unwrap().transform
        = Matrix4D::scaling(0.1, 0.1, 0.1);
    floor.material_mut().reflective = 0.5;

    let mut middle = Shape::sphere();
    *middle.transform_mut() = Matrix4D::translation(-0.5, 1.0, -5.0);
    *middle.material_mut() = Default::default();
    middle.material_mut().color = Color::rgb(1.0, 0.4666, 0.2666);
    middle.material_mut().diffuse = 0.7;
    middle.material_mut().specular = 0.3;
    middle.material_mut().transparency = 0.5;
    middle.material_mut().reflective = 0.5;
    // middle.material.refractive_index = 1.5;

    /*
    middle.material.pattern = Some(Pattern::checker(
        Color::blue(),
        Color::green()
    ));
    */

    let mut right = Shape::capped_cylinder(0.0, 3.0);
    *right.transform_mut() = Matrix4D::translation(1.5, 2.5, -5.5)
        * Matrix4D::scaling(0.25, 0.25, 0.25)
        // * Matrix4D::rotation_x(std::f64::consts::PI / 1.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_y(std::f64::consts::PI / 8.0);
    *right.material_mut() = Default::default();
    // right.material.color = Color::rgb(1.0, 0.6666, 0.2666);
    right.material_mut().diffuse = 0.7;
    right.material_mut().specular = 0.3;
    right.material_mut().reflective = 0.3;
    right.material_mut().transparency = 0.9;
    right.material_mut().refractive_index = 1.5;

    let mut left = Shape::cube();
    *left.transform_mut() = Matrix4D::translation(-1.5, 1.0, -0.75)
        * Matrix4D::scaling(0.33, 0.33, 0.33)
        * Matrix4D::rotation_x(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_y(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 4.0);
    *left.material_mut() = Default::default();
    // left.material.color = Color::rgb(0.8666, 0.2, 0.2);
    left.material_mut().diffuse = 0.7;
    left.material_mut().specular = 0.3;
    left.material_mut().transparency = 0.5;
    left.material_mut().refractive_index = 1.1;

    let mut world = World::empty();
    world.light_source = PointLight::new(
        Color::rgb(1.0, 1.0, 1.0),
        Tuple4D::point(-10.0, 10.0, -10.0),
    );
    world.objects = vec![
        floor,
        middle,
        right,
        left,
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
