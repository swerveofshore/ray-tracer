use std::fs;
use std::ffi::OsStr;
use std::path::Path;

use clap::{ app_from_crate, arg };

use ray_tracer_challenge::scene::{ Scene, SceneJson };
use ray_tracer_challenge::tuple::Tuple4D;
use ray_tracer_challenge::matrix::Matrix4D;
use ray_tracer_challenge::color::Color;
use ray_tracer_challenge::light::PointLight;
use ray_tracer_challenge::shape::Shape;
use ray_tracer_challenge::world::World;
use ray_tracer_challenge::camera::Camera;
use ray_tracer_challenge::parallel::parallel_render;
use ray_tracer_challenge::consts::{ CANVAS_WIDTH, CANVAS_HEIGHT };

fn main() {
    let matches = app_from_crate!()
        .arg(
            arg!(
                -s --scene <FILE> "Sets a scene to be rendered."
            )
            .required(false)
            .allow_invalid_utf8(true),
        )
        .arg(
            arg!(
                -o --output <FILE> "Sets the output file (PPM format)."
            )
            .required(false)
            .allow_invalid_utf8(true),
        )
        .arg(
            arg!(
                -j --jobs <VALUE> "Sets the number of threads to be used."
            )
            .required(false),
        )
        .get_matches();

    // Extract the number of jobs, if provided.
    let jobs = match matches.value_of("jobs") {
        Some(j) => j.parse().expect(
            "-j or --jobs argument must be a positive integer."
        ),
        None => 1,
    };

    // Extract the output path, if provided.
    let out = Path::new(
        match matches.value_of_os("output") {
            Some(o) => o,
            None => OsStr::new(ray_tracer_challenge::consts::OUT_FILE),
        }
    );

    // Render from a JSON scene or use a default path.
    if let Some(raw_scene) = matches.value_of_os("scene") {
        let scene_path = Path::new(raw_scene);
        let scene_str = fs::read_to_string(scene_path).expect(
            &format!("Failed to read scene JSON file {:?}.", scene_path)
        );

        let scene_json: SceneJson
            = serde_json::from_str(&scene_str).expect(
                &format!("Failed to parse scene JSON file {:?}", scene_path)
            );

        let scene: Scene = scene_json.into();
        parallel_render(scene.world, scene.camera, jobs, out);
    } else {
        let sphere = Shape::sphere();
        let mut floor = Shape::plane();
        floor.set_transform(Matrix4D::translation(0.0, -4.0, 0.0));

        let mut world = World::empty();
        world.light_source = PointLight::new(
            Color::rgb(0.85, 0.8, 0.65),
            Tuple4D::point(-10.0, 10.0, -10.0),
        );

        world.objects = vec![
            sphere,
            floor,
        ];

        let mut camera = Camera::new(CANVAS_WIDTH, CANVAS_HEIGHT,
            std::f64::consts::PI / 3.0, Matrix4D::identity());
        camera.transform = Matrix4D::view_transform(
            Tuple4D::point(0.0, 1.5, -5.0),
            Tuple4D::point(0.0, 1.0,  0.0),
            Tuple4D::vector(0.0, 1.0, 0.0),
        );

        // Render the world.
        parallel_render(world, camera, jobs, out);
    }
}
