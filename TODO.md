# Future Enhancements

Here are some possible enhancements for this ray tracer.

    * Bounds calculations are repeated frequently, even though they could be
      saved. How could bounds be saved without large modifications to the
      intersection logic?

    * Multithreading should be easy, but it should involve a thread pool of
      sorts (no barriers).

        * The Rust Book discusses a thread pool--perhaps this could be looked
          at.

        * Would some paradigm work better than a straightforward thread pool?

    * Patterns are somewhat boring.

        * Implement a "perturb" function which distorts existing patterns based
          on Perlin (or alternative) noise.

        * Implement a "texture" pattern which applies some image to a solid.

            * Might need to look at image libraries for this.

        * Maybe look at PoVRay for inspiration?

    * More shapes could be implemented.

        * Torus

        * Quadric surfaces

    * Glass and metal defaults could be nice. There's a "glassy sphere" default
      already, and that seems to work nicely. More defaults for metal or "dry"
      materials would be good.

    * Some code is overtly repetitive (for example, `ShapeNode` constructors).

        * Macros may help with condensing this, but should only be used as a
          last-resort.

    * Several functions need better documentation; some notes could be provided
    on math-heavy intersection or lighting functions.

    * Structures should have tests for default values as a sanity check.

    * Test names should be more consistent.
