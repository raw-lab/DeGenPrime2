//! WebAssembly entry points (plain C ABI, no wasm-bindgen needed).
//!
//! JavaScript copies a request into memory obtained from `dgp_alloc`, calls
//! `dgp_call(endpoint, body)`, reads the reply (a little-endian u32 length
//! followed by UTF-8 JSON) and releases every buffer with `dgp_free`.

use std::alloc::{alloc, dealloc, Layout};

#[link(wasm_import_module = "env")]
extern "C" {
    /// Supplied by JavaScript: receives the panic message before the trap.
    fn dgp_panic(ptr: *const u8, len: usize);
}

fn install_panic_hook() {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        std::panic::set_hook(Box::new(|info| {
            let msg = info.to_string();
            unsafe { dgp_panic(msg.as_ptr(), msg.len()) }
        }))
    });
}

#[no_mangle]
pub extern "C" fn dgp_alloc(len: usize) -> *mut u8 {
    unsafe { alloc(Layout::from_size_align_unchecked(len.max(1), 1)) }
}

/// # Safety
/// `ptr` must come from `dgp_alloc(len)` (or a `dgp_call` reply of `len` bytes).
#[no_mangle]
pub unsafe extern "C" fn dgp_free(ptr: *mut u8, len: usize) {
    if !ptr.is_null() {
        dealloc(ptr, Layout::from_size_align_unchecked(len.max(1), 1));
    }
}

/// # Safety
/// Both pointers must reference live buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn dgp_call(ep_ptr: *const u8, ep_len: usize, body_ptr: *const u8, body_len: usize) -> *mut u8 {
    install_panic_hook();
    let ep = std::str::from_utf8(std::slice::from_raw_parts(ep_ptr, ep_len)).unwrap_or("");
    let body = String::from_utf8_lossy(std::slice::from_raw_parts(body_ptr, body_len));
    let out = crate::api::handle(ep, &body);
    let bytes = out.as_bytes();
    let p = dgp_alloc(bytes.len() + 4);
    std::ptr::copy_nonoverlapping((bytes.len() as u32).to_le_bytes().as_ptr(), p, 4);
    std::ptr::copy_nonoverlapping(bytes.as_ptr(), p.add(4), bytes.len());
    p
}
