//! Counts spelled out in prose, against the thing they count.
//!
//! `friction_counts.rs` checks `FRICTION.md`'s own summary sentence against its headings, and that count
//! has not been wrong since. What it does not check is the **six other places** the same total is
//! restated, and a `prose-auditor` pass before the 0.14.0 release found four of them saying twenty-three
//! when the file said twenty-four. Three more documents said fourteen or eighteen scenes when there were
//! twenty-one.
//!
//! Six of that audit's twenty findings were this one shape: a number of scenes or findings written into a
//! document nothing compares against anything. That is a *class*, and a correction closes one instance of
//! it while a test closes the class.
//!
//! # Exact phrases, and why the first design was worse
//!
//! The first version scanned for any number-word near a subject word — "N scenes", "N findings" — and
//! compared every one against the total. It fired three times immediately and **all three were its own
//! fault**: `scenes/README.md` legitimately says "the first of three scenes with nothing to draw", and
//! `AGENTS.md` legitimately says "twelve worked problems" about fifteen example files, because three of
//! them are a quickstart, a benchmark and a README checker rather than worked problems. Subset counts and
//! differently-defined counts are normal prose and a blunt scan cannot tell them from staleness.
//!
//! I had also predicted the wrong failure mode for that design — written down as "the risk is a false
//! pass" — and what it produced was three false failures.
//!
//! So each claim is registered as a **template** with the number left as `{}`, and the test asserts two
//! things: the template filled with the right word is present, and filled with any *other* number-word it
//! is not. That has no false positives, and rewording a sentence makes the test fail with "not found"
//! rather than pass silently — which is the property that matters, because a vacuous pass is how this
//! class of staleness survives in the first place.
//! # Not under `wasm32`
//!
//! Every test here reads a file out of the repository, and a `wasm32-wasip1` runner has no
//! repository — no preopened directory, and nothing at the paths these walk to. They used to live in
//! `pantometry-world`, which the wasm jobs excluded with a flag; the crate moved to `app/` and these
//! two came here instead, where nothing excluded them and five `test` jobs went red.
//!
//! A `cfg` rather than a flag on the job: the reason is a property of the *test* — it needs a
//! filesystem — and a flag on a CI line is a fact about the tests kept somewhere the tests are not.
#![cfg(not(target_family = "wasm"))]

use std::path::{Path, PathBuf};

/// Words for the numbers these documents use. Prose does not say "21".
const WORDS: [&str; 41] = [
    "zero",
    "one",
    "two",
    "three",
    "four",
    "five",
    "six",
    "seven",
    "eight",
    "nine",
    "ten",
    "eleven",
    "twelve",
    "thirteen",
    "fourteen",
    "fifteen",
    "sixteen",
    "seventeen",
    "eighteen",
    "nineteen",
    "twenty",
    "twenty-one",
    "twenty-two",
    "twenty-three",
    "twenty-four",
    "twenty-five",
    "twenty-six",
    "twenty-seven",
    "twenty-eight",
    "twenty-nine",
    "thirty",
    "thirty-one",
    "thirty-two",
    "thirty-three",
    "thirty-four",
    "thirty-five",
    "thirty-six",
    "thirty-seven",
    "thirty-eight",
    "thirty-nine",
    "forty",
];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("crates/pantometry has two ancestors")
        .to_path_buf()
}

/// Assert that `template` filled with `want` appears in `relative`, and that no other number-word fills
/// it.
///
/// A missing file is skipped: this crate is `publish = false` but its tests run from a checkout, and a
/// packaging arrangement that hid a document is not this test's business. A missing *phrase* is a
/// failure, which is the whole point.
fn phrase(relative: &str, template: &str, want: usize) {
    let Ok(text) = std::fs::read_to_string(root().join(relative)) else {
        return;
    };
    assert!(
        template.contains("{}"),
        "the template needs a hole for the number"
    );
    // Case-insensitively, because a claim at the start of a sentence is capitalised and that is a
    // property of English rather than of the count. Two of these phrases are, and the first run of this
    // test failed on both.
    let text = text.to_lowercase();
    let filled = |word: &str| template.to_lowercase().replace("{}", word);
    let wanted = filled(WORDS[want]);
    assert!(
        text.contains(&wanted),
        "{relative} no longer contains {wanted:?}. If the sentence was reworded, update the \
         template here — a phrase that cannot be found is the failure this test is shaped to give \
         instead of a silent pass"
    );
    for (n, word) in WORDS.iter().enumerate() {
        if n == want {
            continue;
        }
        let wrong = filled(word);
        // A shorter number-word is a *substring* of a longer one: "twenty-one worlds described as data"
        // contains "one worlds described as data", and "twenty-four findings" contains "four findings".
        // So a hit only counts when what precedes it is not part of a word — which is one more false
        // positive this design produced before it stopped producing them.
        let spurious = |at: usize| {
            text[..at]
                .chars()
                .next_back()
                .is_some_and(|c| c == '-' || c.is_alphanumeric())
        };
        let mut from = 0;
        while let Some(at) = text[from..].find(&wrong) {
            let hit = from + at;
            assert!(
                spurious(hit),
                "{relative} contains {wrong:?}; the count is {want}"
            );
            from = hit + 1;
        }
    }
}

/// How many findings `FRICTION.md` holds, and how many are fixed — computed the same way
/// `friction_counts.rs` computes them, so the two tests cannot disagree about the number.
fn friction_totals() -> Option<(usize, usize)> {
    let text = std::fs::read_to_string(root().join("app/pantometry-world/FRICTION.md")).ok()?;
    let findings = text
        .lines()
        .filter(|l| {
            l.strip_prefix("## ")
                .and_then(|r| r.split_once('.'))
                .is_some_and(|(n, _)| !n.is_empty() && n.chars().all(|c| c.is_ascii_digit()))
        })
        .count();
    let fixed = text.lines().filter(|l| l.starts_with("**Fixed")).count();
    (findings > 0 && fixed > 0).then_some((findings, fixed))
}

/// **The findings total is the same in all eight places it is written.**
///
/// Four of these were stale at once before the 0.14.0 release — `CLAUDE.md`, two agent files and
/// `FRICTION.md`'s own closing section — while `friction_counts.rs` passed, because it checks one
/// sentence in one file. One file under test is not enough when eight restate the number.
///
/// The eighth was found by reading rather than by this test: `consumer-advocate.md`'s opening
/// paragraph said **twenty-two findings, seventeen fixed** — a pair from an earlier moment left
/// standing as if it were current, while the same file's closing section carried the right one.
/// Guarded now, which is the only reason to have noticed it twice.
#[test]
fn the_findings_total_agrees_everywhere_it_is_written() {
    let Some((findings, fixed)) = friction_totals() else {
        return;
    };
    println!("  {findings} findings, {fixed} fixed");

    phrase(
        "app/pantometry-world/FRICTION.md",
        "**{} of the thirty-four are fixed**",
        fixed,
    );
    phrase(
        "app/pantometry-world/FRICTION.md",
        "{} findings, and the source has shifted",
        findings,
    );
    phrase(
        "CLAUDE.md",
        "{} findings from using the SDK as a stranger",
        findings,
    );
    // The *total*, with the fixed count spelled beside it. Hard-coding that half made this
    // template refuse a correct README the day a finding was actioned -- twice in one commit,
    // because the sentence appears in two shapes. Both halves come from the count now.
    phrase(
        "README.md",
        &format!("**{{}} findings, {} fixed", WORDS[fixed]),
        findings,
    );
    // The other half is the total minus the fixed count, so the two add up by construction rather
    // than by somebody remembering. Hard-coded, it refused a correct README the day finding 8 was
    // actioned -- and this test is what found that sentence at all, because a grep for
    // "twenty-eight" missed it: the fixed count and the argued-down count are spelled in the same
    // line and only one of them was the number being searched for.
    phrase(
        "README.md",
        &format!(
            "findings, {{}} fixed and {} argued down",
            WORDS[findings - fixed]
        ),
        fixed,
    );
    phrase(
        "CHANGELOG.md",
        "has already found {} places it is awkward",
        findings,
    );
    phrase(
        "CHANGELOG.md",
        "places it is awkward, {} of which have been changed",
        fixed,
    );
    phrase(
        "app/pantometry-world/src/lib.rs",
        "beside this crate. {} of the thirty-four are",
        fixed,
    );
    phrase(
        ".claude/agents/README.md",
        "by a distance: {} findings",
        findings,
    );
    phrase(".claude/agents/README.md", "findings, {} fixed", fixed);
    // "Six of the thirty-four" -- both halves, and the first one was hard-coded here while the
    // second came from the count. So this template refused a correct file the moment a finding was
    // actioned, which is the same shape as the two in README.md above and was found the same way.
    phrase(
        ".claude/agents/consumer-advocate.md",
        &format!(
            "{} of the {{}} findings were recorded",
            WORDS[findings - fixed]
        ),
        findings,
    );
    // The opening paragraph, which restates both halves in a different sentence from the closing
    // one. It read "Twenty-two findings ... seventeen of them fixed" while the file's other end
    // said thirty-four and twenty-nine.
    phrase(
        ".claude/agents/consumer-advocate.md",
        &format!(
            "**{{}}** findings have\ncome out of it, {} of them fixed",
            WORDS[fixed]
        ),
        findings,
    );
}

/// **The number of agents is the number of agent files.**
///
/// Two of the seven were describing a workspace nobody could see. `domain-builder`'s *description*
/// — the line an agent picker shows — said the recipe came from "the six existing domains", and
/// its body said "Five exist", while eleven do. Neither is a count this file could have derived
/// from the agents themselves, so it derives them from `crates/`.
///
/// The count of agents is derived from the directory, because that is the set that can be
/// enumerated: adding a ninth file and forgetting the two tables is the failure this shape exists
/// to give instead of a silent pass.
#[test]
fn the_agent_team_counts_itself_and_the_domains_it_describes() {
    let dir = root().join(".claude/agents");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return;
    };
    let agents = entries
        .filter_map(Result::ok)
        .filter(|e| {
            let n = e.file_name().to_string_lossy().into_owned();
            n.ends_with(".md") && n != "README.md"
        })
        .count();
    if agents == 0 {
        return;
    }
    println!("  {agents} agents");
    phrase(
        ".claude/agents/README.md",
        "{}, each built around work",
        agents,
    );
    phrase("CLAUDE.md", "holds {} reviewers", agents);

    // The domains `domain-builder` claims to have learned from, against the crates that are
    // domains: everything in `crates/` that is not the facade, the units, the kernel, or one of
    // the three layers above physics.
    let Ok(crates) = std::fs::read_dir(root().join("crates")) else {
        return;
    };
    let not_a_domain = [
        "pantometry",
        "pantometry-units",
        "pantometry-core",
        "pantometry-scene",
        "pantometry-view",
        "pantometry-shape",
    ];
    let domains = crates
        .filter_map(Result::ok)
        .filter(|e| e.path().is_dir())
        .filter(|e| {
            let n = e.file_name().to_string_lossy().into_owned();
            n.starts_with("pantometry-") && !not_a_domain.contains(&n.as_str())
        })
        .count();
    println!("  {domains} domains");
    phrase(
        ".claude/agents/domain-builder.md",
        "the recipe the {} existing domains established",
        domains,
    );
    phrase(".claude/agents/domain-builder.md", "**{}** exist", domains);
}

/// **The scene count is the same in all seven places it is written.**
///
/// `scene.rs` reads the directory rather than a list, so adding a scene never breaks anything — which is
/// exactly why the sentences about the count drift. Two agent files said fourteen and `scenes/README.md`
/// said eighteen in one place while saying twenty-one in four others.
///
/// The subset claims in the same file — "eleven of these twenty-one have a domain with no field",
/// "twelve of the twenty-one have geometry to export" — are covered here only for the *total* half of
/// each sentence. Their own numerators are not counted by anything and are left alone rather than
/// guarded badly.
#[test]
fn the_scene_count_agrees_everywhere_it_is_written() {
    let dir = root().join("app/pantometry-world/scenes");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return;
    };
    let scenes = entries
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().is_some_and(|x| x == "json"))
        .count();
    if scenes == 0 {
        return;
    }
    println!("  {scenes} scenes");

    let readme = "app/pantometry-world/scenes/README.md";
    phrase(readme, "{} worlds described as data", scenes);
    phrase(readme, "which is what all {} here", scenes);
    // These two sentences also carry a *numerator* — how many scenes have a fieldless domain, how many
    // export geometry — and the templates deliberately cover only the total. A numerator is a different
    // count with a different source, and folding it in here would make the guard fail every time one
    // moved, which is what it did when scene 22 arrived: the template held "Eleven of these {} scenes" and
    // twelve became correct.
    //
    // Both were verified against the binary rather than reasoned about — 12 scenes print a "not drawn"
    // domain and 13 export a `.gltf` — but nothing counts them continuously, and a guard that demanded a
    // template edit for every such move would be a guard somebody rewrites rather than reads.
    phrase(readme, "of these {} scenes", scenes);
    phrase(readme, "of the {} scenes have geometry", scenes);
    phrase(readme, "runs all {} on every commit", scenes);
    // **A sixth place, found by reading rather than by this test.** `CLAUDE.md` says what lives in
    // `app/` and listed "the twenty-eight scenes" while thirty were there — a live claim about the
    // present, not one of the two sentences in `editor-core` that say "the twenty-eight scenes of
    // the time" and mean it. The five above are all in one README, which is how a count can be
    // guarded in five places and still be wrong in a sixth.
    phrase("CLAUDE.md", "and so do the {} scenes", scenes);
    // And a seventh, written in the same session that found the sixth: the editor's README says
    // where a new scene's duration and frame count come from, which is a claim about how many
    // scenes there are to take them from.
    phrase(
        "app/editor-core/README.md",
        "against the {} scenes on disk",
        scenes,
    );
    // The crate table moved to `ARCHITECTURE.md` when `README.md` was cut back to what this is
    // and how to run it. The sentence is the same sentence; this test found the move by failing,
    // which is what it is for.
    phrase(
        "ARCHITECTURE.md",
        "with {} scenes across all eleven domains",
        scenes,
    );
    phrase(
        ".claude/agents/domain-builder.md",
        "{} ship, all run by CI",
        scenes,
    );
}

/// **The example count in `RELEASING.md`'s own counting command is the number of example files.**
///
/// Not a prose claim — a shell snippet, which is worse, because it looks like a measurement. The command
/// there globbed `crates/pantometry/examples/*.rs` and missed the one example that lives in `pantometry-optics`,
/// so a release following it would have counted fourteen of fifteen and had no way to know.
///
/// **`AGENTS.md`'s "twelve worked problems" is deliberately not checked against this.** Three of the
/// fifteen files are a quickstart, a benchmark and a README checker, and calling them worked problems
/// would be the wrong sentence rather than the wrong number. A test that insisted they agree would be
/// asserting a definition, and it fired on exactly that when this file's first draft tried.
#[test]
fn the_releasing_example_command_counts_every_example() {
    let mut examples = 0;
    let mut stack = vec![root().join("crates")];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.filter_map(Result::ok) {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|x| x == "rs")
                && path.parent().is_some_and(|p| p.ends_with("examples"))
            {
                examples += 1;
            }
        }
    }
    if examples == 0 {
        return;
    }
    println!("  {examples} example files across the workspace");
    let Ok(text) = std::fs::read_to_string(root().join("RELEASING.md")) else {
        return;
    };
    let want = format!("# examples: {examples}");
    assert!(
        text.contains(&want),
        "RELEASING.md's counting command should say {want:?}; there are {examples} example files \
         and the answer belongs beside the command that finds them"
    );
    assert!(
        !text.contains("ls crates/pantometry/examples/*.rs | grep -vc common"),
        "the old glob missed the example in pantometry-optics"
    );
}

/// The **six crates that are not a physics**, so the domain count is a subtraction with a stated
/// list rather than a number somebody remembers.
///
/// Naming them beats deriving them. Every plausible derivation is wrong somewhere: "depends on
/// the kernel" catches `pantometry-scene`, `pantometry-view` and the facade; "has a `Domain`
/// impl" is true of test fixtures. A list fails *loudly* when a seventh non-physics crate
/// arrives, which is the direction a guard should fail in.
const NOT_A_PHYSICS: [&str; 6] = [
    "pantometry",       // the facade
    "pantometry-units", // under everything
    "pantometry-core",  // the kernel, which knows no physics by construction
    "pantometry-shape", // layer 0, input
    "pantometry-scene", // layer 2
    "pantometry-view",  // layer 3
];

/// How many crates there are, and how many of them are a physics, wherever prose says so.
///
/// **Added because all three of these had drifted at once.** At 0.18.0 the citation record and
/// the Zenodo deposition both described "sixteen crates ... ten domain crates" — one short on
/// each, and in the two documents a reader outside this repository is most likely to be handed.
/// `ARCHITECTURE.md`'s layer diagram said ten as well. Nothing here read any of them, and the
/// three that this file already covered were the three that had stayed correct.
#[test]
fn the_crate_and_domain_counts_agree_everywhere_they_are_written() {
    let Ok(entries) = std::fs::read_dir(root().join("crates")) else {
        return;
    };
    let names: Vec<String> = entries
        .filter_map(Result::ok)
        .filter(|e| e.path().is_dir())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    if names.is_empty() {
        return;
    }
    for expected in NOT_A_PHYSICS {
        assert!(
            names.iter().any(|n| n == expected),
            "{expected} is in NOT_A_PHYSICS and not in crates/ -- the list has gone stale"
        );
    }
    let crates = names.len();
    let domains = crates - NOT_A_PHYSICS.len();
    println!("  {crates} crates, {domains} of them a physics");

    // The two documents that leave this repository: one is what a citation manager reads, the
    // other is what the DOI record shows.
    for f in ["CITATION.cff", ".zenodo.json"] {
        phrase(f, "A Rust workspace of {} crates", crates);
        phrase(f, "{} domain crates built on it", domains);
    }
    phrase("ARCHITECTURE.md", "+ {} domain crates", domains);
    // **The front page's own two.** `README.md` was cut back to what this is and how to run it,
    // and what survived includes the crate count twice — in the `cargo add` comment and in the
    // line pointing at the map. A short document's numbers are read more, not less.
    phrase("README.md", "all {} published crates", crates);
    phrase("README.md", "the map: three layers, the {} crates", crates);
    phrase(
        "ARCHITECTURE.md",
        "The physics layer is {} crates deep",
        domains,
    );
}
