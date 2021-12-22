use std::rc::Rc;

use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;
// use ray_tracer_challenge::extra::hexagon;
use ray_tracer_challenge::obj::ObjParser;

const CANVAS_WIDTH: usize = 240 * 16;
const CANVAS_HEIGHT: usize = 135 * 16; 
const OBJ_FILE: &'static str = "./models/teapot.obj";
const OUT_FILE: &'static str = "./out.ppm";

fn main() {
    println!("Parsing OBJ file {}...", OBJ_FILE);
    let mut obj_parser = ObjParser::new(OBJ_FILE);
    obj_parser.parse();
    println!("...done.");

    let mut world = World::empty();

    world.light_source = PointLight::new(
        Color::rgb(0.6, 0.8, 1.0),
        Tuple4D::point(-10.0, 10.0, -10.0),
    );

    world.objects = vec![
        Rc::clone(obj_parser.groups.get("").unwrap())
    ];

    let mut camera = Camera::new(CANVAS_WIDTH, CANVAS_HEIGHT,
        std::f64::consts::PI / 3.0, Matrix4D::identity());
    camera.transform = Matrix4D::view_transform(
        Tuple4D::point(0.0, 1.5, -5.0),
        Tuple4D::point(0.0, 1.0,  0.0),
        Tuple4D::vector(0.0, 1.0, 0.0),
    ) * Matrix4D::translation(0.0, -2.0, 8.0);

    println!("Rendering...");
    let canvas = camera.render(&mut world);
    canvas.save(OUT_FILE).unwrap();
    println!("...done.");
    println!("PPM file saved to {}.", OUT_FILE);
}
