use std::rc::Rc;
use std::cell::RefCell;

use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::pattern::Pattern;
use ray_tracer_challenge::shape::Shape;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;
use ray_tracer_challenge::obj::ObjParser;

const CANVAS_WIDTH: usize = 125; // * 8; // 240 * 8;
const CANVAS_HEIGHT: usize = 100; // * 8; // 135 * 8; 
const OBJ_FILE: &'static str = "./models/old-teapot.obj";
const OUT_FILE: &'static str = "./out.ppm";

fn main() {
    println!("Parsing OBJ file {}...", OBJ_FILE);
    let mut obj_parser = ObjParser::new(OBJ_FILE);
    obj_parser.parse();
    println!("...done.");

    // let model = obj_parser.groups.get("").unwrap();
    let mut models: Vec<_> =
        obj_parser.groups.values().map(|g| Rc::clone(&g)).collect();

    /*
    model.borrow_mut().transform
        = Matrix4D::rotation_y(3.0 * std::f64::consts::PI / 4.0);

    for child in model.borrow().children().unwrap() {
        child.borrow_mut().material.transparency = 0.8;
        child.borrow_mut().material.reflective = 0.1;
        child.borrow_mut().material.refractive_index = 1.5;
    }
    */

    let floor = Rc::new(RefCell::new(Shape::plane()));
    floor.borrow_mut().transform
        = Matrix4D::translation(0.0, -4.0, 0.0);
    floor.borrow_mut().material.pattern
        = Some(Pattern::checker(Color::black(), Color::white()));
    floor.borrow_mut().material.reflective = 0.8;

    let cube = Rc::new(RefCell::new(Shape::cube()));
    cube.borrow_mut().transform = Matrix4D::scaling(1.5, 1.5, 1.5)
        * Matrix4D::rotation_y(std::f64::consts::PI / 8.0)
        * Matrix4D::rotation_z(std::f64::consts::PI / 8.0);
    cube.borrow_mut().material.color = Color::rgb(1.0, 1.0, 0.0);

    let sphere = Rc::new(RefCell::new(Shape::sphere()));
    sphere.borrow_mut().transform =
        Matrix4D::rotation_y(std::f64::consts::PI / 8.0)
        * Matrix4D::translation(0.5, 0.0, 0.0);
    sphere.borrow_mut().material.color = Color::red();

    let difference = Shape::csg_difference(cube, sphere);
    difference.borrow_mut().transform = Matrix4D::translation(-5.0, 0.0, 0.0);

    let mut world = World::empty();
    world.light_source = PointLight::new(
        Color::rgb(0.85, 0.8, 0.65),
        Tuple4D::point(-10.0, 10.0, -10.0),
    );

    world.objects = vec![
        Rc::clone(&floor),
        Rc::clone(&difference),
    ];
    world.objects.append(&mut models);

    let mut camera = Camera::new(CANVAS_WIDTH, CANVAS_HEIGHT,
        std::f64::consts::PI / 3.0, Matrix4D::identity());
    camera.transform = Matrix4D::view_transform(
        Tuple4D::point(0.0, 1.5, -5.0),
        Tuple4D::point(0.0, 1.0,  0.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
    ) * Matrix4D::translation(2.0, 0.0, 12.0);

    println!("Rendering...");
    let canvas = camera.render(&mut world);
    canvas.save(OUT_FILE).unwrap();
    println!("...done.");
    println!("PPM file saved to {}.", OUT_FILE);
}
