use serde::{ Serialize, Deserialize };

use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::light::{ PointLight, Material };
use crate::shape::Shape;
use crate::world::World;
use crate::camera::Camera;

#[derive(Debug)]
pub struct Scene {
    pub world: World,
    pub camera: Camera,
}

impl From<SceneJson> for Scene {
    fn from(scene_json: SceneJson) -> Scene {
        // Create the camera transform from the view parameters.
        let camera_transform = Matrix4D::view_transform(
            Tuple4D { w: 1.0, ..(&scene_json.camera_from).into() },
            Tuple4D { w: 1.0, ..(&scene_json.camera_to).into() },
            (&scene_json.camera_up).into()
        );

        // Create the camera.
        let camera = Camera::new(
            scene_json.canvas_width,
            scene_json.canvas_height,
            scene_json.field_of_view,
            camera_transform
        );

        // Create the world.
        let mut world = World::empty();
        world.light_source = scene_json.light.into();
        world.objects = scene_json.shapes.into_iter().map(|x| x.into())
            .collect();

        Scene { world, camera, }
    }
}

#[derive(Serialize, Deserialize)]
pub struct SceneJson {
    canvas_width: usize,
    canvas_height: usize,
    field_of_view: f64,

    camera_from: Vec<f64>,
    camera_to: Vec<f64>,
    camera_up: Vec<f64>,

    light: LightJson,
    shapes: Vec<ShapeJson>,
}

#[derive(Clone, Serialize, Deserialize)]
struct LightJson {
    intensity: Vec<f64>,
    position: Vec<f64>,
}

impl From<LightJson> for PointLight {
    fn from(light_json: LightJson) -> PointLight {
        PointLight {
            intensity: (&light_json.intensity).into(),
            position: Tuple4D { w: 1.0, ..(&light_json.position).into() }
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct ShapeJson {
    ty: String,
    transform: Vec<f64>,
    material: Option<MaterialJson>,
    children: Option<Vec<ShapeJson>>,
}

impl From<ShapeJson> for Shape {
    fn from(shape_json: ShapeJson) -> Shape {
        let mut shape = match shape_json.ty.as_str() {
            // Primitives
            "empty"            => Shape::empty(),
            "sphere"           => Shape::sphere(),
            "plane"            => Shape::plane(),
            "cube"             => Shape::cube(),
            "cylinder"         => Shape::cylinder(),
            "bounded_cylinder" => Shape::bounded_cylinder(-1.0, 1.0),
            "capped_cylinder"  => Shape::capped_cylinder(-1.0, 1.0),
            "bounded_cone"     => Shape::bounded_cone(0.0, 1.0),
            "bounded_dn_cone"  => Shape::bounded_cone(-1.0, 1.0),
            "capped_cone"      => Shape::capped_cone(0.0, 1.0),
            "capped_dn_cone"   => Shape::capped_cone(-1.0, 1.0),

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
        };

        // Set the shape transform
        let mut shape_transform: [f64; 16] = [0.0; 16];
        for i in 0..shape_json.transform.len() {
            shape_transform[i] = shape_json.transform[i];
        }
        shape.set_transform(shape_transform.into());

        // Set the shape material
        if let Some(m) = shape_json.material {
            shape.material = m.into();
        }

        shape
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct MaterialJson {
    color: Vec<f64>,

    // TODO handle patterns
    // pattern: Option<PatternJson>,

    ambient: Option<f64>,
    diffuse: Option<f64>,
    specular: Option<f64>,
    shininess: Option<f64>,

    reflective: Option<f64>,
    refractive_index: Option<f64>,
    transparency: Option<f64>,
}

impl From<MaterialJson> for Material {
    fn from(material_json: MaterialJson) -> Material {
        let mut material = Material {
            color: (&material_json.color).into(),
            ..Default::default()
        };

        if let Some(a) = material_json.ambient {
            material.ambient = a;
        }

        if let Some(d) = material_json.diffuse {
            material.diffuse = d;
        }

        if let Some(s) = material_json.specular {
            material.specular = s;
        }

        if let Some(s) = material_json.shininess {
            material.shininess = s;
        }

        if let Some(r) = material_json.reflective {
            material.reflective = r;
        }

        if let Some(ri) = material_json.refractive_index {
            material.refractive_index = ri;
        }

        if let Some(t) = material_json.transparency {
            material.transparency = t;
        }

        material
    }
}
