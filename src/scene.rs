use serde::{ Serialize, Deserialize };

use crate::world::World;
use crate::camera::Camera;
use crate::shape::Shape;

pub struct Scene {
    pub world: World,
    pub camera: Camera,
}

impl From<SceneJson> for Scene {
    fn from(scene_json: SceneJson) -> Scene {
        // Take the first sixteen elements of the scene_json transform vector,
        // use them to create a Matrix4D for the camera.
        let mut camera_transform: [f64; 16] = Default::default();
        for (i, elem)
            in scene_json.camera_transform.into_iter().take(16).enumerate() {
            camera_transform[i] = elem; 
        }

        // Create the camera.
        let camera = Camera::new(
            scene_json.canvas_width,
            scene_json.canvas_height,
            scene_json.field_of_view,
            camera_transform.into()
        );

        // Create the world.
        let mut world = World::empty();
        world.objects = scene_json.shapes.into_iter().map(|x| x.into()).collect();

        Scene { world, camera, }
    }
}

#[derive(Serialize, Deserialize)]
pub struct SceneJson {
    canvas_width: usize,
    canvas_height: usize,
    field_of_view: f64,
    camera_transform: Vec<f64>,

    light: LightJson,
    shapes: Vec<ShapeJson>,
}

#[derive(Clone, Serialize, Deserialize)]
struct LightJson {
    intensity: Vec<f64>,
    position: Vec<f64>,
}

#[derive(Clone, Serialize, Deserialize)]
struct ShapeJson {
    ty: String,
    transform: Vec<f64>,
    children: Option<Vec<ShapeJson>>,
}

impl From<ShapeJson> for Shape {
    fn from(shape_json: ShapeJson) -> Shape {
        match shape_json.ty.as_str() {
            // Primitives
            "empty" => Shape::empty(),
            "sphere" => Shape::sphere(),
            "plane" => Shape::plane(),
            "cube" => Shape::cube(),
            "cylinder" => Shape::cylinder(),
            "bounded_cylinder" => Shape::bounded_cylinder(-1.0, 1.0),
            "capped_cylinder" => Shape::capped_cylinder(-1.0, 1.0),
            "bounded_cone" => Shape::bounded_cone(0.0, 1.0),
            "bounded_dn_cone" => Shape::bounded_cone(-1.0, 1.0),
            "capped_cone" => Shape::capped_cone(0.0, 1.0),
            "capped_dn_cone" => Shape::capped_cone(-1.0, 1.0),

            // Group-likes
            "group" => {
                let mut group = Shape::group();
                if let Some(children) = shape_json.children {
                    // Convert all the child shape JSONs to shapes.
                    *group.children_mut().unwrap()
                        = children.into_iter().map(|x| x.into()).collect();
                }

                // It's okay to have an empty group (no children).
                group
            },
            "union" => {
                let mut csg_union = Shape::empty();
                if let Some(children) = shape_json.children {
                    if children.len() < 2 {
                        panic!("CSG union must have at least two operands.");
                    }

                    let left = children[0].clone().into();
                    let right = children[1].clone().into();
                    csg_union = Shape::csg_union(left, right);
                }

                csg_union
            },
            "intersection" => {
                let mut csg_intersection = Shape::empty();
                if let Some(children) = shape_json.children {
                    if children.len() < 2 {
                        panic!("CSG union must have at least two operands.");
                    }

                    let left = children[0].clone().into();
                    let right = children[1].clone().into();
                    csg_intersection = Shape::csg_intersection(left, right);
                }

                csg_intersection
            }
            "difference" => {
                let mut csg_difference = Shape::empty();
                if let Some(children) = shape_json.children {
                    if children.len() < 2 {
                        panic!("CSG union must have at least two operands.");
                    }

                    let left = children[0].clone().into();
                    let right = children[1].clone().into();
                    csg_difference = Shape::csg_difference(left, right);
                }

                csg_difference
            },

            // Models
            // TODO
            _ => panic!("Unrecognized shape type in scene description JSON."),
        }
    }
}
