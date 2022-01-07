use std::thread;
use std::sync::mpsc;
use std::sync::{ Arc, Mutex };
use std::path::Path;

use crate::world::World;
use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::consts::REFLECTION_RECURSION_DEPTH;

/// A message sent between threads.
///
/// A message can either be a pixel to be rendered in a `Scene`, or a signal
/// to terminate the completed rendering thread.
pub enum Message {
    Pixel(usize, usize),
    Terminate,
}

/// A worker which holds a rendering thread.
///
/// Each worker runs until a `Terminate` message is sent. Basically, a worker
/// waits for pixels to render, and once it receives a pixel, the worker renders
/// it.
struct Worker {
    /// A thread for rendering pixels.
    thread: Option<thread::JoinHandle<()>>,
}

impl Worker {
    /// Creates a new worker which immediately begins rendering.
    ///
    /// Once a worker begins, it does not stop until it receives a `Terminate`
    /// message. Calling this function allocates a thread which awaits messages
    /// until the thread terminates.
    fn new(world: Arc<World>, camera: Arc<Camera>, canvas: Arc<Mutex<Canvas>>,
        receiver: Arc<Mutex<mpsc::Receiver<Message>>>) -> Worker {

        let thread = thread::spawn(move || loop { 
            // Obtain the message being executed.
            let message: Message = receiver.lock().unwrap().recv().unwrap();
            
            match message {
                Message::Pixel(x, y) => {
                    // Render a pixel on the canvas.
                    let ray = camera.ray_for_pixel(x, y);
                    let color = world.color_at(ray, REFLECTION_RECURSION_DEPTH);
                    canvas.lock().unwrap().write_pixel(x, y, &color);
                },

                Message::Terminate => {
                    // Exit the worker thread loop, terminating the thread.
                    break;
                }
            }
        });

        Worker { thread: Some(thread) }
    }
}

/// A thread pool for rendering a `Scene`.
///
/// Effectively, `n` threads (specified by the user) are allocated up-front as
/// the program executes. Pixels are sent to each worker thread until the entire
/// `Scene` is rendered.
///
/// When a worker finishes rendering a pixel, it awaits another pixel location
/// to render; these locations are sent by the `sender`.
///
/// See function `parallel_render` for the top-level execution logic regarding
/// the renderer.
pub struct ThreadPool {
    /// Workers allocated for rendering.
    workers: Vec<Worker>,

    /// A sender handle for sending messages to each `Worker`.
    sender: mpsc::Sender<Message>,
}

impl ThreadPool {
    /// Creates a new thread pool of size `size` for rendering a scene.
    ///
    /// Each possible pixel on the `canvas` is sent via an `mpsc` channel to
    /// render on each worker thread. The "scene" is a combination of `world`
    /// and `camera`, where rays are cast from the camera onto the world to
    /// produce colors for each pixel.
    ///
    /// These pixels are sent using method `execute`. This function only
    /// instantiates the threads needed for rendering.
    pub fn new(size: usize, world: World, camera: Camera,
        canvas: Arc<Mutex<Canvas>>) -> ThreadPool {
        // There should be at least one thread to run workers.
        assert!(size > 0);

        let (sender, receiver) = mpsc::channel();

        let world = Arc::new(world);
        let camera = Arc::new(camera);
        let receiver = Arc::new(Mutex::new(receiver));

        let mut workers = Vec::with_capacity(size);

        for _ in 0..size {
            workers.push(Worker::new(
                Arc::clone(&world),
                Arc::clone(&camera),
                Arc::clone(&canvas),
                Arc::clone(&receiver)
            ));
        }

        ThreadPool { workers, sender }
    }

    /// Send a pixel location to any listening worker to render said pixel.
    pub fn execute(&mut self, message: Message) {
        self.sender.send(message).unwrap();
    }
}

/// Custom drop implementation for `ThreadPool`.
///
/// When no more pixel locations are being sent to the `ThreadPool`, terminate
/// all threads started by the `ThreadPool`.
///
/// This prevents the thread pool from waiting to join indefinitely.
impl Drop for ThreadPool {
    fn drop(&mut self) {
        // Send the terminate signal to each worker.
        for _ in &self.workers {
            self.sender.send(Message::Terminate).unwrap();
        }

        // Join all of the worker threads.
        for worker in &mut self.workers {
            if let Some(thread) = worker.thread.take() {
                thread.join().unwrap();
            }
        }
    }
}

/// Render a scene in parallel.
///
/// This is the primary driver code for the executable. Using a scene defined
/// by `world` and `camera`, rays are cast by the `camera` onto the `world` to
/// calculate the color of a pixel at each pixel location.
///
/// Pixel colors are stored in a `Canvas`, which is updated across all worker
/// threads in an instantiated `ThreadPool`.
///
/// Once all threads finish rendering all pixel locations, the resultant render
/// is written to file `out` via the `Canvas`.
///
/// # Examples
///
/// Render the default world with a camera distanced from the world origin:
///
/// ```
/// # use ray_tracer_challenge::tuple::Tuple4D;
/// # use ray_tracer_challenge::matrix::Matrix4D;
/// # use ray_tracer_challenge::camera::Camera;
/// # use ray_tracer_challenge::world::World;
/// # use ray_tracer_challenge::parallel::parallel_render;
/// use std::path::Path;
/// use ray_tracer_challenge::consts::{ CANVAS_WIDTH, CANVAS_HEIGHT };
/// 
/// // Create the default world
/// let world: World = Default::default();
/// let view = Matrix4D::view_transform(
///     Tuple4D::point(0.0, 1.5, -5.0),
///     Tuple4D::point(0.0, 1.0, 0.0),
///     Tuple4D::vector(0.0, 1.0, 0.0)
/// );
///
/// // Create a camera for creating a canvas with dimensions 50 by 40
/// let camera = Camera::new(50, 40, std::f64::consts::PI / 3.0, view);
///
/// // Execute a parallel render with two threads, save to file 'out.ppm'
/// parallel_render(world, camera, 2, &Path::new("out.ppm"));
/// ```
pub fn parallel_render(world: World, camera: Camera, jobs: usize, out: &Path) {
    let vsize = camera.vsize;
    let hsize = camera.hsize;
    let canvas = Arc::new(Mutex::new(Canvas::new(hsize, vsize)));

    println!("Rendering using {} threads...", jobs);
    {
        let mut thread_pool = ThreadPool::new(
            jobs, world, camera, Arc::clone(&canvas)
        );

        for y in 0..vsize {
            for x in 0..hsize {
                thread_pool.execute(Message::Pixel(x, y));
            }
        }
    }

    canvas.lock().unwrap().save(out).unwrap();
    println!("...done.\n");
    println!("Saved render to {:?}.", out);
}
