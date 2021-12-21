use std::rc::Rc;
use std::cell::RefCell;

use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::shape::ShapeNode;
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::pattern::Pattern;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;

const CANVAS_WIDTH: usize = 960 * 4;
const CANVAS_HEIGHT: usize = 540 * 4; 

fn main() {
    let mut backdrop = ShapeNode::plane();
    backdrop.transform = Matrix4D::translation(0.0, 0.0, 10.0)
        * Matrix4D::rotation_x(std::f64::consts::PI / 2.0);
    backdrop.material.pattern = Some(
        Pattern::checker(
            Color::black(),
            Color::white()
        )
    );
    backdrop.material.reflective = 0.8;

    let mut parent = ShapeNode::cube();
    parent.transform = Matrix4D::scaling(4.0, 4.0, 4.0)
        * Matrix4D::rotation_x(std::f64::consts::PI / 8.0)
        * Matrix4D::rotation_y(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 8.0);
    parent.material.color = Color::rgb(0.2, 0.2, 0.4);
    parent.material.transparency = 0.5;
    parent.material.reflective = 0.2;
    parent.material.refractive_index = 2.0;

    let mut child_01 = ShapeNode::cube();
    child_01.transform = Matrix4D::scaling(2.0, 2.0, 2.0)
        * Matrix4D::rotation_x(std::f64::consts::PI / 8.0)
        * Matrix4D::rotation_y(std::f64::consts::PI / 4.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 8.0)
        * Matrix4D::translation(3.0, 0.0, 0.0);
    child_01.material.color = Color::rgb(0.4, 0.4, 0.8);
    child_01.material.transparency = 0.5;
    child_01.material.reflective = 0.2;
    child_01.material.refractive_index = 2.0;

    let mut world = World::empty();

    world.light_source = PointLight::new(
        Color::rgb(1.0, 0.8, 0.6),
        Tuple4D::point(-10.0, 10.0, -10.0),
    );

    world.objects = vec![
        Rc::new(RefCell::new(backdrop)),
        Rc::new(RefCell::new(parent)),
        Rc::new(RefCell::new(child_01))
    ];

    let mut camera = Camera::new(CANVAS_WIDTH, CANVAS_HEIGHT,
        std::f64::consts::PI / 3.0, Matrix4D::identity());
    camera.transform = Matrix4D::view_transform(
        Tuple4D::point(0.0, 1.5, -5.0),
        Tuple4D::point(0.0, 1.0,  0.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
    ) * Matrix4D::translation(0.0, 0.0, 48.0);

    let canvas = camera.render(&mut world);
    canvas.save("out.ppm").unwrap(); }
