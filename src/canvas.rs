use std::io;
use std::io::Write;
use std::fs::File;
use std::path::Path;

use crate::color::Color;

/// A canvas for drawing pixels.
///
/// This structure mostly stores the results of the ray tracer. Once the user
/// specifies the desired image width and height, the `Camera` generates rays
/// which are cast onto a `World`'s `Object`s.
///
/// The canvas stores the resulant colors for each pixel ray. Once execution
/// finishes, the `Canvas` can be used to save the pixels to an image file.
///
/// For now, only PPM images are supported.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct Canvas {
    /// The width of the canvas, in pixels.
    pub width: usize,

    /// The height of the canvas, in pixels.
    pub height: usize,

    /// The pixels of the canvas, stored as a flattened vector.
    pixels: Vec<Color>,
}

impl Canvas {
    /// Creates a new canvas with specified width and height.
    ///
    /// This function allocates a `Vec<Color>` of size `width * height`, which
    /// may take up a decent amount of memory, depending on image size.
    pub fn new(width: usize, height: usize) -> Canvas {
        Canvas {
            width,
            height,
            pixels: vec![Color::rgb(0.0, 0.0, 0.0); width * height]
        }
    }

    /// Saves a canvas to a PPM file.
    ///
    /// Lines in the PPM file are clamped to 70 columns. If some color exceeds
    /// the 70 column mark on a line, it is moved to the next line over.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        // Create a new file
        let mut out = File::create(path)?;

        // Write PPM header to file, as well as metadata
        writeln!(&mut out, "P3")?;
        writeln!(&mut out, "{} {}", self.width, self.height)?;
        writeln!(&mut out, "255")?; // Maximum color value

        // Write pixels to file, making sure that no line exceeds 70 columns
        let mut col = 1;
        for pixel in self.pixels.iter() {
            // Scale pixels to 255.0, then convert to strings
            let r_mod = (pixel.r * 255.0).clamp(0.0, 255.0).ceil() as usize;
            let g_mod = (pixel.g * 255.0).clamp(0.0, 255.0).ceil() as usize;
            let b_mod = (pixel.b * 255.0).clamp(0.0, 255.0).ceil() as usize;
            let r_str = r_mod.to_string();
            let g_str = g_mod.to_string();
            let b_str = b_mod.to_string();

            // Check if any color surpasses the 70 column marker
            if col + r_str.len() > 70 {
                write!(&mut out, "\n{} {} {}", r_str, g_str, b_str)?;
                col = r_str.len() + g_str.len() + b_str.len() + 3;
            } else if col + r_str.len() + g_str.len() > 70 {
                write!(&mut out, " {}\n{} {}", r_str, g_str, b_str)?;
                col = g_str.len() + b_str.len() + 2;
            } else if col + r_str.len() + g_str.len() + b_str.len() > 70 {
                write!(&mut out, " {} {}\n{}", r_str, g_str, b_str)?;
                col = b_str.len() + 1;
            // Otherwise, write colors as normal
            } else {
                if col != 1 {
                    write!(&mut out, " ")?;
                    col += 1;
                }

                write!(&mut out, "{} {} {}", r_str, g_str, b_str)?;
                col += r_str.len() + g_str.len() + b_str.len() + 2;
            }
        }

        // Terminate the PPM file with a newline
        write!(&mut out, "\n")?;

        // If no errors occur, return Ok
        Ok(())
    }

    /// Writes a color to a location on the `Canvas`.
    ///
    /// Out-of-bounds pixels are ignored. Pixels are specified in row-column
    /// order, where `y` is the row of the pixel, and `x` is the column. Rows
    /// and columns are zero-indexed.
    ///
    /// # Examples
    ///
    /// Writing a pixel to the fourth column, second row on an 8-by-8 canvas:
    ///
    /// ```
    /// # use ray_tracer_challenge::color::Color;
    /// # use ray_tracer_challenge::canvas::Canvas;
    /// let purple = Color::rgb(1.0, 0.0, 1.0);
    /// let mut canvas = Canvas::new(8, 8);
    /// canvas.write_pixel(4, 2, &purple);
    /// assert_eq!(canvas.read_pixel(4, 2).unwrap(), purple);
    /// ```
    pub fn write_pixel(&mut self, x: usize, y: usize, pixel: &Color) {
        // Silently ignore out-of-bounds pixels
        if x >= self.width || y >= self.height {
            return;
        }

        self.pixels[(y * self.width) + x] = *pixel;
    }

    /// Reads a color from a location on the `Canvas`.
    ///
    /// Pixels are specified in row-column order, where `y` is the row of the
    /// pixel, and `x` is the column. Rows and columns are zero-indexed. If
    /// the specified pixel location is out-of-bounds, `None` is returned by
    /// this function.
    ///
    /// See method `write_pixel` for an example of writing and reading a pixel
    /// from a `Canvas`.
    pub fn read_pixel(&self, x: usize, y: usize) -> Option<Color> {
        // Return nothing if pixel is out-of-bounds
        if x >= self.width || y >= self.height {
            return None
        }

        Some(self.pixels[(y * self.width) + x])
    }
}
