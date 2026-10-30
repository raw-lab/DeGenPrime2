//! Gives `degenprime.exe` its icon and version details on Windows.
//!
//! The resources are prebuilt by `scripts/make_windows_resource.sh`, so a
//! build needs no resource compiler and no build dependencies. Builds for
//! every other platform are unaffected.

use std::path::Path;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=assets/windows");
    let var = |k: &str| std::env::var(k).unwrap_or_default();
    if var("CARGO_CFG_TARGET_OS") != "windows" || std::env::var_os("CARGO_FEATURE_CLI").is_none() {
        return;
    }
    let file = match (var("CARGO_CFG_TARGET_ENV").as_str(), var("CARGO_CFG_TARGET_ARCH").as_str()) {
        // link.exe and lld-link take compiled resources directly, for any architecture.
        ("msvc", _) => "degenprime.res",
        // The MinGW linker takes them as a COFF object.
        ("gnu", "x86_64") => "degenprime-x86_64-gnu.o",
        _ => return,
    };
    let path = Path::new(&var("CARGO_MANIFEST_DIR")).join("assets").join("windows").join(file);
    if path.is_file() {
        println!("cargo:rustc-link-arg-bin=degenprime={}", path.display());
    } else {
        println!("cargo:warning={} is missing, so degenprime.exe will have no icon", path.display());
    }
}
