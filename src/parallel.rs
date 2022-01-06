use std::thread;
use std::sync::mpsc;
use std::sync::{ Arc, Mutex };

use crate::world::World;
use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::consts::{ NUM_THREADS, OUT_FILE, REFLECTION_RECURSION_DEPTH };

pub enum Message {
    Pixel(usize, usize),
    Terminate,
}

struct Worker {
    thread: Option<thread::JoinHandle<()>>,
}

impl Worker {
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

pub struct ThreadPool {
    workers: Vec<Worker>,
    sender: mpsc::Sender<Message>,
}

impl ThreadPool {
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

    pub fn execute(&mut self, message: Message) {
        self.sender.send(message).unwrap();
    }
}

impl Drop for ThreadPool {
    fn drop(&mut self) {
        for _ in &self.workers {
            self.sender.send(Message::Terminate).unwrap();
        }

        for worker in &mut self.workers {
            if let Some(thread) = worker.thread.take() {
                thread.join().unwrap();
            }
        }
    }
}

pub fn parallel_render(world: World, camera: Camera) {
    let vsize = camera.vsize;
    let hsize = camera.hsize;
    let canvas = Arc::new(Mutex::new(Canvas::new(hsize, vsize)));

    println!("Rendering using {} threads...", NUM_THREADS);
    {
        let mut thread_pool = ThreadPool::new(
            NUM_THREADS, world, camera, Arc::clone(&canvas)
        );

        for y in 0..vsize {
            for x in 0..hsize {
                thread_pool.execute(Message::Pixel(x, y));
            }
        }
    }

    canvas.lock().unwrap().save(OUT_FILE).unwrap();
    println!("...done.\n");
    println!("Saved render to {}.", OUT_FILE);
}
