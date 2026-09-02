//! **A heading that moves takes its links with it, and nothing said so.**
//!
//! `README.md` was 883 lines and most of it was explanation. Cutting it back to what this is and
//! how to run it moved nine sections into `EVIDENCE.md`, two into `ARCHITECTURE.md` and one into
//! `CONTRIBUTING.md` — and one sentence in the moved text said *"see the section above"* with a
//! same-file anchor. Above, in a file it was no longer in. Markdown does not complain, GitHub
//! renders the link, and clicking it does nothing.
//!
//! So: every `](#anchor)` in the documents that ship, against the headings of the file it is in,
//! and every `](FILE.md#anchor)` against the headings of that file.
//!
//! # What this cannot see
//!
//! A link to a file that exists but is the wrong one, and a link to a *web* address. The first
//! needs a reader; the second needs the network, which no test here is allowed to want.
//!
//! # Not under `wasm32`
//!
//! Every document is read off a disk, and a `wasm32` target has none: with nothing readable the
//! walk finds no links and the "a check that stopped finding links would pass forever" assertion
//! fires, which is the assertion doing its job about the wrong thing. `counts_in_prose` and
//! `citation_is_valid` carry the same line for the same reason, and this one shipped without it —
//! the local gate has no wasm step and CI's `test (wasm32-wasip1, wasmtime)` job found it.

#![cfg(not(target_family = "wasm"))]

/// The repository root, from this crate's manifest.
fn root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("crates/pantometry has two ancestors")
        .to_path_buf()
}

/// GitHub's rule, near enough: lowercase, punctuation dropped, spaces to hyphens.
fn slug(heading: &str) -> String {
    let h = heading.trim_start_matches('#').trim().to_lowercase();
    let kept: String = h
        .chars()
        .filter(|c| {
            !matches!(
                c,
                '`' | '*' | '[' | ']' | '(' | ')' | ',' | '.' | ':' | '\'' | '"'
            )
        })
        .collect();
    let mut out = String::new();
    let mut hyphen = false;
    for c in kept.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c);
            hyphen = false;
        } else if !hyphen && !out.is_empty() {
            out.push('-');
            hyphen = true;
        }
    }
    out.trim_end_matches('-').to_string()
}

/// The `#`-headings of a file, as anchors.
fn anchors(text: &str) -> std::collections::BTreeSet<String> {
    text.lines()
        .filter(|l| l.starts_with('#'))
        .map(slug)
        .collect()
}

/// Every `](...)` link in `text`, as `(file, anchor)` with an empty file meaning "this one".
fn links(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let bytes: Vec<char> = text.chars().collect();
    let mut i = 0;
    while i + 1 < bytes.len() {
        if bytes[i] == ']' && bytes[i + 1] == '(' {
            let mut j = i + 2;
            let mut target = String::new();
            while j < bytes.len() && bytes[j] != ')' && bytes[j] != ' ' {
                target.push(bytes[j]);
                j += 1;
            }
            if let Some((file, anchor)) = target.split_once('#') {
                if !anchor.is_empty() && !target.starts_with("http") {
                    out.push((file.to_string(), anchor.to_string()));
                }
            }
            i = j;
        }
        i += 1;
    }
    out
}

/// The documents a reader arrives at, and the two beside the things they describe.
const DOCS: [&str; 12] = [
    "README.md",
    "AGENTS.md",
    "ARCHITECTURE.md",
    "EVIDENCE.md",
    "EXAMPLES.md",
    "CONTRIBUTING.md",
    "CHANGELOG.md",
    "CLAUDE.md",
    "RELEASING.md",
    "app/README.md",
    "app/editor-core/README.md",
    "app/pantometry-world/FRICTION.md",
];

/// **Every anchor a document links to is a heading that exists.**
#[test]
fn every_anchor_a_document_links_to_is_a_heading_somewhere() {
    let root = root();
    let mut checked = 0;
    let mut broken = Vec::new();
    for doc in DOCS {
        let Ok(text) = std::fs::read_to_string(root.join(doc)) else {
            continue;
        };
        let here = anchors(&text);
        for (file, anchor) in links(&text) {
            checked += 1;
            let (label, set) = if file.is_empty() {
                (doc.to_string(), here.clone())
            } else {
                let at = root.join(doc).parent().map_or_else(
                    || root.join(&file),
                    |dir| {
                        let beside = dir.join(&file);
                        if beside.is_file() {
                            beside
                        } else {
                            root.join(&file)
                        }
                    },
                );
                match std::fs::read_to_string(&at) {
                    Ok(t) => (file.clone(), anchors(&t)),
                    Err(_) => {
                        broken.push(format!("{doc}: {file} does not exist"));
                        continue;
                    }
                }
            };
            if !set.contains(&anchor) {
                broken.push(format!("{doc}: {label}#{anchor} is not a heading there"));
            }
        }
    }
    // A check that stopped finding links would pass forever, so the walk has to say it walked.
    // **Three, measured** — the documents cross-reference each other by *file* far more than by
    // heading, and only three links name an anchor at all. That is a thin surface and worth
    // saying: this catches the anchor that moved, and it would not have caught nine of the ten
    // sections moving if none of them had been linked to.
    assert!(
        checked > 0,
        "no anchored links found in any document — either the walk broke or they all went"
    );
    println!("  {checked} anchored links");
    assert!(broken.is_empty(), "{}", broken.join("\n"));
}
