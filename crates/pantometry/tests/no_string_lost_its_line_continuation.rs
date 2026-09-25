//! **No string in the tree has lost its line continuation.**
//!
//! A long Rust string is written across lines with a backslash, which takes the newline and the
//! next line's indentation with it:
//!
//! ```text
//! "... a scene holding one would \
//!  not load"                          is   "... a scene holding one would not load"
//! ```
//!
//! Edited through a Python string, the backslash and the newline are Python's continuation and go
//! first, so the indentation stays inside the Rust string and the two source lines become one:
//! `would` and `not` twenty spaces apart, on a line longer than rustfmt's hundred columns. rustfmt
//! does not reflow a string literal, so nothing said so. **Thirty-nine lines** had it, from
//! twenty-nine commits going back to 2026-08-08 — nine of them in `src/`, including the refusal a
//! person gets for asking for a GPU and the bytes `pantometry-view` writes for a surface panel.
//!
//! # The rule, and why it is not wider
//!
//! A string literal holding a run of nine or more spaces between two words, **on a line over a
//! hundred columns**. The length is what separates the defect from a table: fourteen lines in the
//! tree hold such a run on purpose — `fit`'s header, the examples' "name … value" columns — and
//! every one is short, because it was written as one line. A lost continuation is two lines joined,
//! each within the width, so the join is not.
#![cfg(not(target_family = "wasm"))]

use std::path::{Path, PathBuf};

/// Every `.rs` file under the directories a person edits, skipping build output.
fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let name = entry.file_name();
        if path.is_dir() {
            if name != "target" && name != ".git" && name != "node_modules" {
                rust_files(&path, out);
            }
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

/// The contents of every double-quoted string on a line, escapes respected.
///
/// A character literal holding a quote, `'"'`, is stepped over so it does not open a string.
fn strings(line: &str) -> Vec<String> {
    let chars: Vec<char> = line.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == '\'' && i + 2 < chars.len() && chars[i + 1] == '"' && chars[i + 2] == '\'' {
            i += 3;
            continue;
        }
        if chars[i] != '"' {
            i += 1;
            continue;
        }
        let mut body = String::new();
        i += 1;
        while i < chars.len() && chars[i] != '"' {
            if chars[i] == '\\' && i + 1 < chars.len() {
                body.push(chars[i]);
                body.push(chars[i + 1]);
                i += 2;
            } else {
                body.push(chars[i]);
                i += 1;
            }
        }
        out.push(body);
        i += 1;
    }
    out
}

/// Whether `s` holds nine or more spaces with something other than a space on both sides.
fn has_run(s: &str) -> bool {
    let chars: Vec<char> = s.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == ' ' {
            let start = i;
            while i < chars.len() && chars[i] == ' ' {
                i += 1;
            }
            if i - start >= 9 && start > 0 && i < chars.len() {
                return true;
            }
        } else {
            i += 1;
        }
    }
    false
}

#[test]
fn no_string_lost_its_line_continuation() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut files = Vec::new();
    for top in ["crates", "app", "bindings", "tools"] {
        rust_files(&root.join(top), &mut files);
    }
    // **Two hundred and fifty-three are tracked**, so a walk that found a handful went somewhere
    // else, and a check over nothing passes in exactly the way this one must not.
    assert!(
        files.len() > 200,
        "found {} Rust files under {}; the walk is not looking at this repository",
        files.len(),
        root.display()
    );

    let mut found = Vec::new();
    for path in &files {
        let Ok(text) = std::fs::read_to_string(path) else {
            continue;
        };
        for (n, line) in text.lines().enumerate() {
            if line.chars().count() <= 100 || line.trim_start().starts_with("//") {
                continue;
            }
            if strings(line).iter().any(|s| has_run(s)) {
                let shown = path.strip_prefix(&root).unwrap_or(path);
                found.push(format!("{}:{}", shown.display(), n + 1));
            }
        }
    }
    println!("  {} Rust files read", files.len());
    assert!(
        found.is_empty(),
        "{} line(s) hold a string with a run of spaces on a line over 100 columns — a line \
         continuation lost to an edit, which reads as words far apart:\n  {}",
        found.len(),
        found.join("\n  ")
    );
}

/// The two helpers find what they are for, on the kind of line that was fixed, and pass a table.
#[test]
fn the_rule_sees_a_join_and_passes_a_column() {
    let joined = format!("\"a scene holding one would{}not load\"", " ".repeat(14));
    assert!(strings(&joined).iter().any(|s| has_run(s)));
    let column = "\"drawn               {:>10} runs\"";
    assert!(
        strings(column).iter().any(|s| has_run(s)),
        "a column has the run too"
    );
    assert!(
        column.chars().count() <= 100,
        "and it is the line's length, not the run, that lets it through"
    );
    let quote_char = "let c = '\"'; let s = \"no run here\";";
    assert_eq!(strings(quote_char), vec!["no run here".to_string()]);
}
