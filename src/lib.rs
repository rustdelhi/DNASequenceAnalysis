pub mod aliner;
pub mod mutation_detection;
pub mod reader;

/// Initiate tracing
pub fn init_logging() -> tracing_appender::non_blocking::WorkerGuard {
    let (non_blocking, guard) = tracing_appender::non_blocking(std::io::stdout());
    tracing_subscriber::fmt().with_writer(non_blocking).init();
    guard
}
