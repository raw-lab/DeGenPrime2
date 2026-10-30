=== DeGenPrime ===
Contributors: rawlab
Tags: bioinformatics, pcr, primer design, dna, science
Requires at least: 5.8
Requires PHP: 7.4
Stable tag: 2.3.0
License: Non-commercial (see LICENSE)

Degenerate PCR primer design from DNA sequences, aligned or not, running in the visitor's browser.

== Description ==

DeGenPrime designs PCR primers that work across a group of related sequences. Paste an
alignment (FASTA or CLUSTAL), or sequences that are not aligned yet (MAFFT is built in), and it returns ranked primer pairs with IUPAC degenerate codes,
melting temperatures, binding sites in every sequence, a simulated gel and CSV/FASTA exports.
It can also check primers you already have.

The engine is the same Rust code as the DeGenPrime command-line tool, compiled to WebAssembly.
It runs inside the visitor's browser: sequences are never sent to your server, and the plugin
adds no database tables, options or cron jobs.

== Installation ==

1. Plugins > Add New > Upload Plugin, choose degenprime-wordpress-plugin.zip, then Install and Activate.
2. Tools > DeGenPrime > "Create the Primer design page", or add [degenprime] to any page.

Shortcode attributes: view (design, results, check, manual, about) and layout (contained, wide, full).

== Frequently Asked Questions ==

= Are sequences uploaded to my site? =
No. All calculations happen in the visitor's browser.

= The app says the engine could not start =
A Content-Security-Policy header must allow 'wasm-unsafe-eval' (script-src) and blob: (worker-src),
and optimisation plugins should not combine or delay degenprime.js. See Tools > DeGenPrime.

= Can it be listed in the WordPress.org plugin directory? =
Not under its current licence: the directory requires GPL-compatible plugins and DeGenPrime is
licensed for non-commercial use. Install it from the zip instead.

== Changelog ==

= 2.3.0 =
* Built from DeGenPrime 2.3.0: unaligned sequences are aligned in the browser with the built-in
  pure-Rust MAFFT, larger inputs, a Cancel button and the DeGenPrime logo.

= 2.2.0 =
* First WordPress release, built from DeGenPrime 2.2.0.
