use std::io;
use std::io::Write;
use std::fs::File;
use std::path::Path;

use crate::color::Color;

#[derive(Clone, Default, Debug, PartialEq)]
pub struct Canvas {
    pub width: usize,
    pub height: usize,

    pixels: Vec<Color>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Canvas {
        Canvas {
            width,
            height,
            pixels: vec![Color::rgb(0.0, 0.0, 0.0); width * height]
        }
    }

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

    pub fn write_pixel(&mut self, x: usize, y: usize, pixel: &Color) {
        // Silently ignore out-of-bounds pixels
        if x >= self.width || y >= self.height {
            return;
        }

        self.pixels[(y * self.width) + x] = *pixel;
    }

    pub fn read_pixel(&self, x: usize, y: usize) -> Option<Color> {
        // Return nothing if pixel is out-of-bounds
        if x >= self.width || y >= self.height {
            return None
        }

        Some(self.pixels[(y * self.width) + x])
    }
}
