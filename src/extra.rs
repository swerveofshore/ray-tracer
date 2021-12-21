#![allow(unused)]

use std::rc::Rc;
use std::cell::RefCell;

use crate::tuple::Tuple4D;
use crate::matrix::Matrix4D;
use crate::shape::{ ShapeNode, ShapeType };

#[derive(Copy, Clone, Default)]
pub struct Projectile {
    pub pos: Tuple4D,
    pub vel: Tuple4D,
}

#[derive(Copy, Clone, Default)]
pub struct Environment {
    pub grav: Tuple4D,
    pub wind: Tuple4D,
}

impl Environment {
    /// Advances a projectile by one tick, based on environmental conditions.
    ///
    /// The projectile returned from this function is the parameter `proj`
    /// subjected to one "tick" of time in the environment `env`.
    /// 
    /// Velocity, gravity and wind are all represented as vectors.
    pub fn tick(&self, proj: &Projectile) -> Projectile {
        let pos = proj.pos + proj.vel;
        let vel = proj.vel + self.grav + self.wind;

        Projectile { pos, vel }
    }
}

pub fn hexagon_corner() -> ShapeNode {
    let mut corner = ShapeNode::sphere();
    corner.transform = Matrix4D::translation(0.0, 0.0, -1.0)
        * Matrix4D::scaling(0.25, 0.25, 0.25);

    corner
}

pub fn hexagon_edge() -> ShapeNode {
    let mut edge = ShapeNode::cylinder();
    if let ShapeType::Cylinder(ref mut min, ref mut max, _) = edge.ty {
        *min = 0.0;
        *max = 1.0;
    }

    edge.transform = Matrix4D::translation(0.0, 0.0, -1.0)
        * Matrix4D::rotation_y(-std::f64::consts::PI / 6.0)
        * Matrix4D::rotation_z(-std::f64::consts::PI / 2.0)
        * Matrix4D::scaling(0.25, 1.0, 0.25);

    edge
}

/*
// TODO: How to insert a ShapeNode in a World, when the World needs to own the
// ShapeNode? (Can't use Rc<T>s without making the World use them).
pub fn hexagon() -> ShapeNode {
    let hex = Rc::new(RefCell::new(ShapeNode::group()));
}
*/
