/// Tuples and arithmetic on tuples.
///
/// Structure `Tuple4D` is used frequently throughout this crate to represent
/// vectors or points in space. Most important methods on vectors and points
/// are included in this module.
pub mod tuple;

/// Matrices and arithmetic on matrices.
///
/// More specifically, 4x4 matrices. Most if not all of the ray tracer logic
/// can be represented with 4x4 matrices, so only a `Matrix4D` is exposed to
/// the user.
///
/// Some other matrices (such as 2x2 and 3x3 matrices) are defined in this
/// module, however they are private. This may change in a later release.
pub mod matrix;

/// Rays and arithmetic on rays.
///
/// As implied by the definition of a `Ray4D`, a ray is defined by an origin
/// (a point in space) and a direction (a vector pointing "out" from the
/// origin).
///
/// Mostly, rays are used for intersecting `Shape`s. See the `shape`,
/// `intersect` and `geometry` modules for more information.
pub mod ray;

/// Shapes, shape types and logic regarding shapes.
///
/// Most of these shapes are provided as a convenience; shapes such as
/// `ShapeType::Sphere` and `ShapeType::Cube` are included for situations where
/// a model is insufficient. (Models are loaded through the `obj` module, which
/// produces either `ShapeType::Triangle`s or `ShapeType::SmoothTriangle`s).
///
/// Shapes can either be primitive (like spheres, cubes or triangles), or
/// group-like (explicitly, `ShapeType::Group`, loosely, constructive solid
/// geometry (CSG) shapes). Primitives hold a transform and their type only.
/// Group-like shapes can hold child shapes, which use the transform of their
/// parent group.
///
/// Several methods are provided for finding the point of intersection on
/// shapes (given a ray), or the normal vector at a particular point on the
/// shapes. These methods are all used for producing renders.
pub mod shape;

/// Intersections and details associated with intersections.
///
/// Defines structures which carry the `t` of intersection (assuming some
/// intersection ray), as well as which object was intersected.
///
/// Other structures are provided for recording intersection computations.
pub mod intersect;

/// Geometry records and bounding boxes.
///
/// Defines structures for associating information with certain geometric
/// primitives (mostly triangles).
///
/// Also contains most bounding box logic.
pub mod geometry;

/// Lights and materials.
///
/// Defines lighting and material records, as well as how to calculate lighting
/// at a given point for other given environmental parameters.
pub mod light;

/// Colors and operations on colors.
///
/// Mostly just a highly specific ordered triple for operations on colors.
pub mod color;

/// Patterns to be applied to objects.
///
/// Similar to module `light`, this module helps define what color is produced
/// when a ray intersects points on a shape.
pub mod pattern;

/// Worlds which rays are traced in.
///
/// Defines a container for shapes and other environmental parameters for which
/// rays may be cast onto a "scene" (not to be confused with crate `scene`,
/// which provides logic for defining scenes).
///
/// Ultimately, a `World` is a light source and a set of objects; the camera and
/// canvas help decide which pixels are drawn.
pub mod world;

/// Cameras for generating a canvas.
///
/// The camera helps define *which* pixels, or which points, are preserved in
/// frame.
///
/// As an aside, the `Camera` struct is *able* to render scenes with method
/// `render`, but this isn't desirable; look at the `parallel` module's function
/// `parallel_render` instead.
pub mod camera;

/// Canvases for generating output images.
///
/// Handles converting the pixels rendered by the camera to a image viewable
/// by an image viewer. As of now, only PPM images can be generated.
pub mod canvas;

/// Scenes and scene descriptions.
///
/// A "scene" is a `World` and a `Camera` to spectate the world.
///
/// This module mostly defines records which can be serialized and deserialized
/// into/from JSON. These records are kept separate from internal structures
/// such as e.g. `Material` to prevent oddities in serialization, but this may
/// change later.
///
/// Either way, this module handles most of the logic regarding turning a user's
/// scene description JSON file into an actual scene which can be rendered.
pub mod scene;

/// OBJ parsing and loading.
///
/// Takes `.obj` files and turns them into `Shape`s, usually composed of
/// multiple groups with several triangles.
pub mod obj;

/// Extra functions and logic from the Ray Tracer Challenge book.
///
/// Most items are irrelevant to the actual ray tracer logic, but they can be
/// somewhat fun (e.g. making a hexagon from `Shape` groups and primitives).
///
/// I did not make the logic for these functions; the functions in this module
/// mostly come from the Ray tracer Challenge Book.
pub mod extra;

/// Constants.
///
/// Some constants are used as defaults, such as `OUT_FILE`. Others are used
/// for reference, such as `GLASS_RI`, in case a user or programmer wants to
/// use the refraction index of glass in their scene description.
pub mod consts;

/// Parallel execution for rendering.
///
/// Since each pixel can be rendered independent of one another in a ray tracer,
/// the functions in this module do exactly that: render each pixel independent
/// of each other in a thread pool.
///
/// The user specifies the amount of threads desired, and the logic in this
/// module executes those threads, and renders pixels on demand (whenever one
/// pixel finishes rendering, another can begin immediately).
///
/// Other optimizations are likely possible for future improvements.
pub mod parallel;

/// Checks if two floating point numbers are near each other (equal).
///
/// Due to error in floating point calculations, two numbers which are
/// effectively equal may be determined unequal. This function considers two
/// floating point numbers "equal" if they are within `consts::FEQ_EPSILON` of
/// each other.
pub fn feq(left: f64, right: f64) -> bool {
    (left - right).abs() < consts::FEQ_EPSILON
}
