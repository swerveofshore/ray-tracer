use std::io::{ self, prelude::* };
use std::path::{ Path, PathBuf };
use std::fs::File;
use std::collections::BTreeMap;

use crate::tuple::Tuple4D;
use crate::shape::Shape;

/// A face in an OBJ file.
type ObjFace = Vec<(usize, Option<usize>, Option<usize>)>;

/// A parser for OBJ files.
#[derive(Clone, Debug)]
pub struct ObjParser {
    pub path: PathBuf,
    pub ignored_lines: usize,

    pub vertices: Vec<Tuple4D>,
    pub normals: Vec<Tuple4D>,
    pub faces: Vec<ObjFace>,
    pub groups: BTreeMap<String, Shape>,
}

impl ObjParser {
    /// Creates a new `ObjParser` to parse the OBJ file at `path_str`.
    pub fn new(path_str: &str) -> ObjParser {
        let mut groups = BTreeMap::new();
        groups.insert("".into(), Shape::group());

        ObjParser {
            path: Path::new(path_str).into(),
            ignored_lines: 0,

            vertices: Vec::new(),
            normals: Vec::new(),
            faces: Vec::new(),
            groups,
        }
    }

    /// Parses an OBJ file.
    ///
    /// Unsupported commands are ignored. Each ignored command increments value
    /// `ignored_lines` by 1.
    pub fn parse(&mut self) {
        let obj_file = File::open(&self.path)
            .expect(&format!("Failed to open file {:?}.", self.path));
        let obj_file_lines = io::BufReader::new(obj_file).lines();

        // Maintain the name of the group that vertices/faces are being
        // inserted into.
        let mut current_group = "".into();
        for line in obj_file_lines {
            // Make sure that the line could be read.
            let line = line.expect(
                &format!("Could not read line from file {:?}.", self.path)
            );

            // Ignore empty lines.
            if line.is_empty() {
                continue;
            }

            // Handle the command. If the command is a `g` command, change the
            // current group to its argument.
            self.handle_command(&line, &mut current_group);
        }
    }

    /// Parses a line of an OBJ file.
    ///
    /// A sample OBJ file may look like the following:
    ///
    /// ```obj
    /// v -1 1 0
    /// v -1 0 0
    /// v 1 0 0
    /// v 1 1 0
    ///
    /// g FirstGroup
    /// f 1 2 3
    /// g SecondGroup
    /// f 1 3 4
    /// ```
    ///
    /// The first letter of each line stands for a function, and all elements
    /// after the first letter are arguments to said function.
    ///
    /// Function `v` specifies a **v**ertex in space, as a point. Function `f`
    /// specifies a **f**ace, composed of vertices.
    ///
    /// Function `g` specifies a named **g**roup, which contains several faces.
    /// Note that groups do not nest; in the above example, `SecondGroup` is
    /// separate from `FirstGroup` (i.e. `SecondGroup` is not a child group).
    fn handle_command(&mut self, line: &str, current_group: &mut String) {
        let params: Vec<&str> = line.split(' ').collect();
        if params[0] == "v" {
            // If there aren't enough parameters to specify a vertex, ignore
            // this line. (Vertices are points in 3D space).
            if params.len() < 4 {
                self.ignored_lines += 1;
                return;
            }

            // TODO implement better error handling
            self.vertices.push(
                Tuple4D::point(
                    params[1].parse().unwrap(),
                    params[2].parse().unwrap(), 
                    params[3].parse().unwrap()
                )
            );
        } else if params[0] == "vn" {
            if params.len() < 4 {
                self.ignored_lines += 1;
                return;
            }

            // TODO implement better error handling
            self.normals.push(
                Tuple4D::vector(
                    params[1].parse().unwrap(),
                    params[2].parse().unwrap(),
                    params[3].parse().unwrap()
                )
            );
        } else if params[0] == "f" {
            // Collect the face specified by this line.
            let mut face: ObjFace = Vec::new();
            for vertex in params.into_iter().skip(1) {
                // If the face involves vertex indices only, parse those.
                if let Ok(v) = vertex.parse() {
                    face.push((v, None, None));
                }
                // Otherwise, extract vertex, texture and normal indices.
                else {
                    let attributes: Vec<_> = vertex.split('/').collect();

                    // If there are more than 3 attributes, ignore this line.
                    if attributes.len() > 3 {
                        self.ignored_lines += 1;
                        return;
                    }

                    // TODO implement better error handling (don't unwrap,
                    // raise an error or ignore at truly malformed attributes).
                    face.push((
                        attributes[0].parse().unwrap(),
                        attributes[1].parse().ok(),
                        attributes[2].parse().ok()
                    ));
                }
            }

            // If any vertex index of the face is larger than the amount of
            // registered vertices, same with each normal index, ignore line.
            if face.iter().any(|&vi| { vi.0 > self.vertices.len() || 
                vi.2.unwrap_or(0) > self.normals.len() }) {
                self.ignored_lines += 1;
                return;
            }

            // Build triangles from the face, and add said triangles to the
            // currently-specified group.
            let triangles = self.fan_triangulation(&face);
            for triangle in triangles {
                if let Some(group) = self.groups.get_mut(current_group) {
                    group.add_child(triangle);
                } else {
                    panic!("OBJ group should be instantiated.");
                }
            }
        } else if line.starts_with("g") {
            // If no group name is provided, ignore this line.
            if params.len() < 2 {
                self.ignored_lines += 1;
                return;
            }

            // If a group has already been instantiated with this name, just
            // change the current group for later commands.
            if self.groups.contains_key(params[1]) {
                *current_group = params[1].into();
            }
            // Otherwise, create a new group with the specified name.
            else {
                self.groups.insert(params[1].into(), Shape::group());
                *current_group = params[1].into();
            }
        }

        // If this line has an unrecognized function, ignore it.
        self.ignored_lines += 1;
    }

    /// Partitions a list of vertices into triangles.
    ///
    /// In an OBJ file, faces can be specified like so:
    ///
    /// ```obj
    /// f 1 2 3 4 5
    /// ```
    ///
    /// This face references vertices 1, 2, 3, 4 and 5. Typically, these
    /// vertices would make a pentagon of some sort.
    ///
    /// Since we only have Triangle primitives, this function takes a
    /// complicated shape like a pentagon and converts it into a list of
    /// Triangles.
    ///
    /// This is done with something called a "fan triangulation." Observe the
    /// following diagram:
    ///
    /// ```text
    ///         B *
    ///          / \
    ///         /   \
    ///        /     \
    ///     A *       * C
    ///       |       |
    ///       |       |
    ///       |       |
    ///     E * ----- * D  
    /// ```
    ///
    /// Albeit somewhat ugly, we can create triangles by taking a starting
    /// vertex (e.g. `A`) and finding two vertices *across* from that vertex
    /// repeatedly, until all vertices have been seen.
    ///
    /// For example, if we start at `A`, we can make a triangle `A-B-C`, then
    /// a triangle `A-C-D`, and finally a triangle `A-D-E`. The first vertex
    /// starts the same, and "rides" the edges of the shape, eventually creating
    /// all possible triangles from the planar form.
    ///
    /// The triangles produced via this method are returned from this function.
    fn fan_triangulation(&self, face: &ObjFace) -> Vec<Shape> {
        let mut triangles = Vec::new();
        
        // Note that the book uses one-based indexing for vertices; this
        // implementation uses zero-based indexing.
        for i in 1..(face.len() - 1) {
            match (face[0].2, face[i].2, face[i+1].2) {
                (Some(n1), Some(n2), Some(n3)) => {
                    triangles.push(
                        Shape::smooth_triangle(
                            // The face contains one-based indices into the
                            // vertices.
                            self.vertices[face[0].0   - 1],
                            self.vertices[face[i].0   - 1],
                            self.vertices[face[i+1].0 - 1],

                            // Each normal is an index as well. TODO implement
                            // better error handling for missing normals.
                            self.normals[n1 - 1],
                            self.normals[n2 - 1],
                            self.normals[n3 - 1]
                        )
                    );
                },

                _ => {
                    triangles.push(
                        Shape::triangle(
                            // The face contains one-based indices into the
                            // vertices.
                            self.vertices[face[0].0   - 1],
                            self.vertices[face[i].0   - 1],
                            self.vertices[face[i+1].0 - 1]
                        )
                    );
                }
            }
        }

        triangles
    }
}

#[test]
fn ignoring_unrecognized_lines() {
    let mut obj_parser = ObjParser::new("./models/gibberish.obj");
    obj_parser.parse();

    assert_eq!(obj_parser.ignored_lines, 5);
}

#[test]
fn vertex_records() {
    let mut obj_parser = ObjParser::new("./models/vertices.obj");
    obj_parser.parse();

    assert_eq!(obj_parser.vertices[0], Tuple4D::point(-1.0, 1.0, 0.0));
    assert_eq!(obj_parser.vertices[1], Tuple4D::point(-1.0, 0.5, 0.0));
    assert_eq!(obj_parser.vertices[2], Tuple4D::point( 1.0, 0.0, 0.0));
    assert_eq!(obj_parser.vertices[3], Tuple4D::point( 1.0, 1.0, 0.0));
}

#[test]
fn parsing_triangle_faces() {
    let mut obj_parser = ObjParser::new("./models/vertices-and-faces.obj");
    obj_parser.parse();

    let default_group = obj_parser.groups.get("").unwrap();
    let children = default_group.children().unwrap();

    let t1 = children[0].triangle_info().unwrap();
    let t2 = children[1].triangle_info().unwrap();

    assert_eq!(t1.p1, obj_parser.vertices[0]);
    assert_eq!(t1.p2, obj_parser.vertices[1]);
    assert_eq!(t1.p3, obj_parser.vertices[2]);
    assert_eq!(t2.p1, obj_parser.vertices[0]);
    assert_eq!(t2.p2, obj_parser.vertices[2]);
    assert_eq!(t2.p3, obj_parser.vertices[3]);
}

#[test]
fn triangulating_polygons() {
    let mut obj_parser = ObjParser::new("./models/vertices-and-polygon.obj");
    obj_parser.parse();

    let default_group = obj_parser.groups.get("").unwrap();
    let children = default_group.children().unwrap();

    let t1 = children[0].triangle_info().unwrap();
    let t2 = children[1].triangle_info().unwrap();
    let t3 = children[2].triangle_info().unwrap();

    assert_eq!(t1.p1, obj_parser.vertices[0]);
    assert_eq!(t1.p2, obj_parser.vertices[1]);
    assert_eq!(t1.p3, obj_parser.vertices[2]);
    assert_eq!(t2.p1, obj_parser.vertices[0]);
    assert_eq!(t2.p2, obj_parser.vertices[2]);
    assert_eq!(t2.p3, obj_parser.vertices[3]);
    assert_eq!(t3.p1, obj_parser.vertices[0]);
    assert_eq!(t3.p2, obj_parser.vertices[3]);
    assert_eq!(t3.p3, obj_parser.vertices[4]);
}

#[test]
fn triangles_in_groups() {
    let mut obj_parser = ObjParser::new("./models/triangles-with-groups.obj");
    obj_parser.parse();

    assert!(obj_parser.groups.contains_key("FirstGroup"));
    assert!(obj_parser.groups.contains_key("SecondGroup"));
}

#[test]
fn vertex_normal_records() {
    let mut obj_parser = ObjParser::new("./models/normal-record.obj");
    obj_parser.parse();

    assert_eq!(obj_parser.normals[0], Tuple4D::vector(0.0, 0.0, 1.0));
    assert_eq!(obj_parser.normals[1], Tuple4D::vector(0.707, 0.0, -0.707));
    assert_eq!(obj_parser.normals[2], Tuple4D::vector(1.0, 2.0, 3.0));
}

#[test]
fn faces_with_normals() {
    let mut obj_parser = ObjParser::new("./models/faces-with-normals.obj");
    obj_parser.parse();

    let default_group = obj_parser.groups.get("").unwrap();
    let children = default_group.children().unwrap();

    let t1 = children[0].smooth_triangle_info().unwrap();
    let t2 = children[1].smooth_triangle_info().unwrap();

    assert_eq!(t1.triangle_info.p1, obj_parser.vertices[0]);
    assert_eq!(t1.triangle_info.p2, obj_parser.vertices[1]);
    assert_eq!(t1.triangle_info.p3, obj_parser.vertices[2]);
    assert_eq!(t1.n1, obj_parser.normals[2]);
    assert_eq!(t1.n2, obj_parser.normals[0]);
    assert_eq!(t1.n3, obj_parser.normals[1]);
    assert_eq!(t1, t2);
}
