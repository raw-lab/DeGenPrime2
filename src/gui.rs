//! The desktop app: `degenprime gui` serves the web interface from memory on
//! 127.0.0.1 and opens it in the default browser. Standard library only: no
//! Python, Tk, Qt or WebView to install, so it behaves the same on Linux,
//! macOS and Windows. Requests run on the native engine through the same JSON
//! API the website uses; nothing leaves the computer.
//!
//! Safety rails: the server binds to loopback only, rejects requests whose
//! `Host` header is not this server (DNS-rebinding defence), and requires a
//! per-session token on every API call (so other web pages cannot use it).
//! A fixed set of worker threads serves connections, so a burst of requests
//! queues instead of starting a thread each.

use crate::api;
use std::collections::hash_map::RandomState;
use std::collections::HashMap;
use std::hash::{BuildHasher, Hasher};
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::{mpsc, Arc, Mutex};
use std::time::Duration;

const INDEX_HTML: &str = include_str!("../web/index.html");
const MAX_BODY: usize = 256 * 1024 * 1024;
/// Connections served at once; more wait in the queue.
const WORKERS: usize = 8;

/// Files served besides the page itself: (path, content type, bytes).
const ASSETS: &[(&str, &str, &[u8])] = &[
    ("/assets/degenprime.js", "text/javascript; charset=utf-8", include_bytes!("../web/assets/degenprime.js")),
    ("/assets/degenprime.css", "text/css; charset=utf-8", include_bytes!("../web/assets/degenprime.css")),
    (
        "/assets/fonts/recursive-latin-mono-normal.woff2",
        "font/woff2",
        include_bytes!("../web/assets/fonts/recursive-latin-mono-normal.woff2"),
    ),
    ("/assets/logo.webp", "image/webp", include_bytes!("../web/assets/logo.webp")),
    ("/assets/mark.png", "image/png", include_bytes!("../web/assets/mark.png")),
    ("/assets/icon-192.png", "image/png", include_bytes!("../web/assets/icon-192.png")),
    ("/assets/icon-512.png", "image/png", include_bytes!("../web/assets/icon-512.png")),
    ("/assets/apple-touch-icon.png", "image/png", include_bytes!("../web/assets/apple-touch-icon.png")),
    ("/assets/favicon.ico", "image/x-icon", include_bytes!("../web/assets/favicon.ico")),
    ("/favicon.ico", "image/x-icon", include_bytes!("../web/assets/favicon.ico")),
    ("/manifest.webmanifest", "application/manifest+json", include_bytes!("../web/manifest.webmanifest")),
];

/// Preferred port; a free one is used if it is taken.
pub const DEFAULT_PORT: u16 = 8737;

pub struct GuiOptions {
    pub port: u16,
    pub open_browser: bool,
}

struct State {
    index: String,
    token: String,
    port: u16,
}

/// 128 random bits for the session token, from the OS where we can read it.
fn random_token() -> String {
    let mut bytes = [0u8; 16];
    #[cfg(unix)]
    let have = std::fs::File::open("/dev/urandom").and_then(|mut f| f.read_exact(&mut bytes)).is_ok();
    #[cfg(not(unix))]
    let have = false;
    if !have {
        // RandomState is seeded from the OS's random source on every platform.
        for (i, chunk) in bytes.chunks_mut(8).enumerate() {
            let mut h = RandomState::new().build_hasher();
            h.write_u64(i as u64 ^ u64::from(std::process::id()));
            chunk.copy_from_slice(&h.finish().to_le_bytes());
        }
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// Start the server and block until the process is stopped.
pub fn serve(opts: GuiOptions) -> Result<(), String> {
    let listener = TcpListener::bind(("127.0.0.1", opts.port))
        .or_else(|_| TcpListener::bind(("127.0.0.1", 0)))
        .map_err(|e| format!("could not open a local port: {e}"))?;
    let port = listener.local_addr().map_err(|e| e.to_string())?.port();
    let token = random_token();
    let config = format!("<script>window.DGP_CONFIG={{\"engine\":\"native\",\"token\":\"{token}\"}};</script>");
    let state = Arc::new(State { index: INDEX_HTML.replace("<!--DGP_CONFIG-->", &config), token, port });
    let url = format!("http://127.0.0.1:{port}/");
    println!("{}", crate::banner());
    println!("  DeGenPrime is running at {url}");
    println!("  Everything stays on this computer. Press Ctrl+C to quit.");
    println!("  Unaligned sequences are aligned with the built-in MAFFT (pure Rust).");
    if opts.open_browser && !open_browser(&url) {
        println!("  Open {url} in your browser.");
    }

    let (tx, rx) = mpsc::sync_channel::<TcpStream>(256);
    let rx = Arc::new(Mutex::new(rx));
    for _ in 0..WORKERS {
        let (rx, st) = (Arc::clone(&rx), Arc::clone(&state));
        std::thread::spawn(move || loop {
            let next = match rx.lock() {
                Ok(guard) => guard.recv(),
                Err(_) => return,
            };
            match next {
                Ok(s) => {
                    let _ = handle(s, &st);
                }
                Err(_) => return,
            }
        });
    }
    for stream in listener.incoming().flatten() {
        if tx.send(stream).is_err() {
            break;
        }
    }
    Ok(())
}

fn handle(stream: TcpStream, st: &State) -> std::io::Result<()> {
    stream.set_read_timeout(Some(Duration::from_secs(30)))?;
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Ok(());
    }
    let mut parts = line.split_whitespace();
    let method = parts.next().unwrap_or("").to_string();
    let target = parts.next().unwrap_or("/").to_string();
    let mut headers: HashMap<String, String> = HashMap::new();
    loop {
        let mut h = String::new();
        if reader.read_line(&mut h)? == 0 {
            break;
        }
        let h = h.trim_end();
        if h.is_empty() {
            break;
        }
        if let Some((k, v)) = h.split_once(':') {
            headers.insert(k.trim().to_ascii_lowercase(), v.trim().to_string());
        }
    }
    let mut out = stream;
    let host = headers.get("host").map(String::as_str).unwrap_or("");
    if host != format!("127.0.0.1:{}", st.port) && host != format!("localhost:{}", st.port) {
        return respond(&mut out, 403, "text/plain; charset=utf-8", b"Forbidden host");
    }
    let len: usize = headers.get("content-length").and_then(|v| v.trim().parse().ok()).unwrap_or(0);
    if len > MAX_BODY {
        return respond(&mut out, 413, "text/plain; charset=utf-8", b"Request too large");
    }
    if headers.get("expect").is_some_and(|v| v.eq_ignore_ascii_case("100-continue")) {
        // curl and some clients wait for this before sending a large body.
        out.write_all(b"HTTP/1.1 100 Continue\r\n\r\n")?;
    }
    let mut body = vec![0u8; len];
    reader.read_exact(&mut body)?;
    let path = target.split('?').next().unwrap_or("/");
    match (method.as_str(), path) {
        ("GET", "/") | ("GET", "/index.html") => {
            respond(&mut out, 200, "text/html; charset=utf-8", st.index.as_bytes())
        }
        ("GET", "/assets/example_anfD.clust") => {
            respond(&mut out, 200, "text/plain; charset=utf-8", crate::EXAMPLE_ALIGNMENT.as_bytes())
        }
        ("GET", p) => match ASSETS.iter().find(|a| a.0 == p) {
            Some((_, ctype, bytes)) => respond(&mut out, 200, ctype, bytes),
            None => respond(&mut out, 404, "text/plain; charset=utf-8", b"Not found"),
        },
        ("POST", p) if p.starts_with("/api/") => {
            if headers.get("x-dgp-token").map(String::as_str) != Some(st.token.as_str()) {
                return respond(
                    &mut out,
                    403,
                    "application/json",
                    br#"{"ok":false,"error":"missing or wrong session token"}"#,
                );
            }
            let body = String::from_utf8(body).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
            let json = api::handle(&p[5..], &body);
            respond(&mut out, 200, "application/json; charset=utf-8", json.as_bytes())
        }
        _ => respond(&mut out, 404, "text/plain; charset=utf-8", b"Not found"),
    }
}

fn respond(out: &mut TcpStream, code: u16, ctype: &str, body: &[u8]) -> std::io::Result<()> {
    let reason = match code {
        200 => "OK",
        204 => "No Content",
        403 => "Forbidden",
        404 => "Not Found",
        413 => "Payload Too Large",
        _ => "Error",
    };
    let head = format!(
        "HTTP/1.1 {code} {reason}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\nCache-Control: no-store\r\nX-Content-Type-Options: nosniff\r\nReferrer-Policy: no-referrer\r\nConnection: close\r\n\r\n",
        body.len()
    );
    out.write_all(head.as_bytes())?;
    out.write_all(body)?;
    out.flush()
}

fn open_browser(url: &str) -> bool {
    use std::process::{Command, Stdio};
    #[cfg(target_os = "macos")]
    let mut cmd = {
        let mut c = Command::new("/usr/bin/open");
        c.arg(url);
        c
    };
    #[cfg(target_os = "windows")]
    let mut cmd = {
        let mut c = Command::new("cmd");
        c.args(["/C", "start", "", url]);
        c
    };
    #[cfg(not(any(target_os = "macos", target_os = "windows")))]
    let mut cmd = {
        let mut c = Command::new("xdg-open");
        c.arg(url);
        c
    };
    cmd.stdout(Stdio::null()).stderr(Stdio::null()).spawn().is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_are_random_hex() {
        let (a, b) = (random_token(), random_token());
        assert_eq!(a.len(), 32);
        assert!(a.chars().all(|c| c.is_ascii_hexdigit()));
        assert_ne!(a, b);
    }

    #[test]
    fn every_asset_is_present() {
        for (path, ctype, bytes) in ASSETS {
            assert!(!bytes.is_empty(), "{path} is empty");
            assert!(!ctype.is_empty());
        }
        assert!(INDEX_HTML.contains("<!--DGP_CONFIG-->"));
    }
}
