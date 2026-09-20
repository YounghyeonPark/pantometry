//! The version this crate tells people to depend on, checked against the version it is.
//!
//! `AGENTS.md` and `README.md` both print an install line. Those lines are prose, and prose
//! stating a number that nothing checks is how this repository has shipped stale counts more
//! than once — the FRICTION totals, the example count, the test count, and `AGENTS.md` itself
//! saying the Python bindings were "not yet on PyPI" for a day after they were.
//!
//! A release bumps `Cargo.toml` and has no reason to touch a markdown file, so this is the
//! failure mode with the least friction of all of them: nothing goes wrong, the number is simply
//! a release behind, and the first person to notice is someone who copied it.

/// Where the repository root is, relative to this crate.
///
/// `None` when the files are absent, which is the packaged case: `cargo package` builds the
/// crate from a tarball that contains no `AGENTS.md`, and a test that failed there would make
/// publishing impossible for a reason that has nothing to do with the crate. Absent is
/// therefore skipped and *present but wrong* is a failure — the distinction that matters.
fn repo_file(name: &str) -> Option<String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(name);
    std::fs::read_to_string(path).ok()
}

/// The `major.minor` a caller should write, which is what the install lines quote.
fn series() -> String {
    let v = env!("CARGO_PKG_VERSION");
    let mut parts = v.split('.');
    let major = parts.next().expect("a version has a major");
    let minor = parts.next().expect("a version has a minor");
    format!("{major}.{minor}")
}

/// **Every `pantometry = "x.y"` in the documentation is the version this crate actually is.**
///
/// Matches the dependency line wherever it appears rather than at a fixed location, so moving
/// the snippet does not silently switch the check off — which is the other way a check like this
/// dies.
#[test]
fn the_install_lines_quote_the_current_series() {
    let want = series();
    let needle = "pantometry = \"";
    let mut checked = 0;

    for file in ["AGENTS.md", "README.md", "CONTRIBUTING.md", "CLAUDE.md"] {
        let Some(text) = repo_file(file) else {
            continue;
        };
        for (line_no, line) in text.lines().enumerate() {
            let Some(at) = line.find(needle) else {
                continue;
            };
            let rest = &line[at + needle.len()..];
            let quoted = rest.split('"').next().unwrap_or("");
            // Only the bare series form is on trial. A line pinning an exact patch, or a path
            // dependency with a version beside it, is saying something else deliberately.
            if quoted.chars().filter(|c| *c == '.').count() != 1 {
                continue;
            }
            checked += 1;
            assert_eq!(
                quoted,
                want,
                "{file}:{} says pantometry = {quoted:?} but this crate is {}. Bump the prose \
                 with the release.",
                line_no + 1,
                env!("CARGO_PKG_VERSION")
            );
        }
    }

    // A check that stopped finding anything would pass forever. If the install lines are
    // reworded out of this shape, this should fail and be rewritten rather than quietly retire.
    if repo_file("AGENTS.md").is_some() {
        assert!(
            checked > 0,
            "no `pantometry = \"x.y\"` line found in the documentation — either it was reworded, \
             in which case update this test, or it was deleted, in which case a caller no \
             longer has one to copy"
        );
    }
}

/// **The citation block's version, its version DOI and its concept DOI are all somebody else's
/// numbers.**
///
/// `README.md` renders a BibTeX entry that a reader copies into a bibliography, and every field in
/// it restates something written elsewhere:
///
/// - `version` is the crate's, which a release bumps in nine places and has no reason to touch
///   here;
/// - `doi` is that version's row in `RELEASING.md`'s table, which is the only place the minted
///   numbers are recorded;
/// - `url` is the **concept** DOI, which never moves and is in `CITATION.cff`.
///
/// This was the last claim in the repository standing on nothing. `RELEASING.md` says so in as
/// many words -- the block "shipped stale for the length of one release" and its `version` "was
/// left at 0.19.0 through the 0.20.0 bump because no list covered it" -- and the note in the
/// block itself records the same thing about 0.16.0 and 0.17.0. A document that has been wrong
/// twice and knows it is a document to check rather than to read.
///
/// The `doi` is the one field a release cannot bump with the rest: it does not exist until the
/// tag has fired the webhook. So this asks `RELEASING.md` rather than asking for a rule, and the
/// order of operations stays what it is.
#[test]
fn the_citation_block_quotes_the_version_and_the_dois_it_should() {
    let Some(readme) = repo_file("README.md") else {
        return; // packaged build
    };
    // The field's value, from `name     = {value},`.
    let field = |name: &str| -> Option<String> {
        readme.lines().find_map(|l| {
            let l = l.trim_start();
            let rest = l
                .strip_prefix(name)?
                .trim_start()
                .strip_prefix('=')?
                .trim_start();
            let inner = rest.strip_prefix('{')?;
            Some(inner.split('}').next()?.to_string())
        })
    };

    let stated = field("version").expect("the BibTeX block names a version");
    assert_eq!(
        stated,
        env!("CARGO_PKG_VERSION"),
        "README.md's BibTeX says version = {{{stated}}} and this crate is {}",
        env!("CARGO_PKG_VERSION")
    );

    let Some(releasing) = repo_file("RELEASING.md") else {
        return;
    };
    // **The table, by this version's row.** `| 0.21.0 | `10.5281/zenodo.22760497` -- ...`. Read by
    // the row rather than by position, so adding a release does not shift the answer.
    let row = format!("| {} |", env!("CARGO_PKG_VERSION"));
    let minted = releasing
        .lines()
        .find(|l| l.trim_start().starts_with(&row))
        .and_then(|l| l.split('`').nth(1).map(str::to_string));
    match minted {
        Some(minted) => {
            let doi = field("doi").expect("the BibTeX block names a doi");
            assert_eq!(
                doi,
                minted,
                "README.md's BibTeX cites {doi} for {} and RELEASING.md's table says {minted}",
                env!("CARGO_PKG_VERSION")
            );
        }
        // **Between the bump and the tag there is no row**, and that is a real state rather than a
        // failure: the DOI does not exist until the release webhook has fired. It is not silence
        // either -- the release procedure is where the edit is listed, and this says which half of
        // it has happened.
        None => println!(
            "  no version-DOI row for {} in RELEASING.md yet -- the `doi` is the edit that comes \
             after the tag",
            env!("CARGO_PKG_VERSION")
        ),
    }

    // The concept DOI, which is the one that never moves. `CITATION.cff` is what a citation
    // manager reads and the BibTeX is what a person copies; the two naming different records
    // would be the worst version of this being wrong.
    let concept = repo_file("CITATION.cff").and_then(|cff| {
        cff.lines().find_map(|l| {
            l.trim_start()
                .strip_prefix("doi:")
                .map(|v| v.trim().to_string())
        })
    });
    if let Some(concept) = concept {
        let url = field("url").expect("the BibTeX block names a url");
        assert!(
            url.ends_with(&concept),
            "README.md's BibTeX `url` is {url} and CITATION.cff's concept DOI is {concept}"
        );
    }
}

/// **The subagents' own statements of the current version are current.**
///
/// `invariant-guard` has a section on the one invariant that cannot be fixed after the fact — a
/// published version is permanent — and it opens by stating what is on crates.io and what the
/// tree is. That line was a release behind at 0.3.0 and a release behind again at 0.4.0: the
/// document whose subject is checking things was the last place in the repository stating a
/// version with nothing standing behind it.
///
/// Only the *tree's* version is checked. What is on crates.io is a fact about the outside world
/// that no test here can know, and the agent is told to run `curl` for it rather than trust the
/// prose — which is the right division, and the reason this checks one number and not two.
#[test]
fn the_agents_know_what_version_the_tree_is() {
    let full = env!("CARGO_PKG_VERSION");
    let needle = "the tree is ";
    let mut checked = 0;

    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../.claude/agents");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return; // packaged build: the agents are not part of any crate
    };
    let mut files: Vec<_> = entries.filter_map(|e| e.ok().map(|e| e.path())).collect();
    files.sort(); // deterministic order, as everything here must be

    for path in files {
        let Ok(text) = std::fs::read_to_string(&path) else {
            continue;
        };
        for (line_no, line) in text.lines().enumerate() {
            let Some(at) = line.find(needle) else {
                continue;
            };
            let stated: String = line[at + needle.len()..]
                .chars()
                .take_while(|c| c.is_ascii_digit() || *c == '.')
                .collect();
            let stated = stated.trim_end_matches('.');
            if stated.is_empty() {
                continue;
            }
            checked += 1;
            assert_eq!(
                stated,
                full,
                "{}:{} says the tree is {stated:?}, and it is {full:?}",
                path.file_name().unwrap_or_default().to_string_lossy(),
                line_no + 1
            );
        }
    }

    assert!(
        checked > 0,
        "no agent states what version the tree is — if that sentence was reworded, update this \
         test rather than letting it pass on finding nothing"
    );
}
