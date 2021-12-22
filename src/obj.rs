use std::io::{ self, prelude::* };
use std::fs::File;
use std::rc::Rc;
use std::cell::RefCell;

use crate::shape::{ ShapeNode, ShapePtr };

#[derive(Clone, Debug, PartialEq)]
enum ObjCommand {
    Vertex(f64, f64, f64),
    Face(usize, usize, usize),
    Group(String),
}

#[derive(Clone, Debug)]
pub struct Obj {
    pub file_name: String,
    pub shape: ShapePtr,
    pub ignored_lines: usize,

    commands: Vec<ObjCommand>,
}

impl Obj {
    pub fn new(file_name: String) -> Obj {
        Obj {
            file_name,
            commands: Vec::new(),
            shape: Rc::new(RefCell::new(ShapeNode::group())),
            ignored_lines: 0,
        }
    }

    pub fn parse(&mut self) {
        let obj_file = File::open(&self.file_name)
            .expect(&format!("Failed to open file {}.", self.file_name));
        let obj_file_lines = io::BufReader::new(obj_file).lines();

        for line in obj_file_lines {
            let line = line.expect(
                &format!("Could not read line from file {}.", self.file_name)
            );

            if let Some(command) = Self::parse_line(&line) {
                self.commands.push(command);
            } else {
                self.ignored_lines += 1;
            }
        }
    }

    /// TODO actually parse each OBJ line
    fn parse_line(line: &str) -> Option<ObjCommand> {
        if line.starts_with("v") {
            Some(ObjCommand::Vertex(0.0, 0.0, 0.0))
        } else if line.starts_with("f") {
            Some(ObjCommand::Face(1, 1, 1))
        } else if line.starts_with("g") {
            Some(ObjCommand::Group("Null".into()))
        } else {
            None
        }
    }

    pub fn shape(&mut self) -> ShapePtr {
        Rc::clone(&self.shape)
    }
}
