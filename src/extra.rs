use crate::tuple::Tuple4D;

#[allow(unused)]
#[derive(Copy, Clone, Default)]
pub struct Projectile {
    pub pos: Tuple4D,
    pub vel: Tuple4D,
}

#[allow(unused)]
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
