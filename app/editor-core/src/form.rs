//! The model tree, and the settings form each of its rows opens.
//!
//! The scene's JSON was the only place a person could add a cooled face or a heat source: the
//! inspector could change a number the file already had and could not give the file a new one. A
//! simulation tool is set up from a tree of what the model holds and a form for whichever part of
//! it is selected, so that is what this is — **over the same text**. The text stays the model.
//! Every change here is one of [`crate::edit`]'s byte splices, so the file keeps its formatting and
//! its key order, and every one can be taken back.
//!
//! # Built from the shape of the file, not from a list of domains
//!
//! A row is made for every domain the scene states, placed or not, and under it a row for every
//! list of objects it holds — `cooling`, `regions`, `nodes`, `perfusions` — with one per entry. A
//! scene-level map whose values are objects (`materials`, `poses`) is a group of its own. Nothing
//! here names a domain kind, so a domain written out of tree gets the same tree. The one place the
//! format is consulted by kind is for the entry a list's **add** starts from, which is
//! [`pantometry_world::templates::ENTRIES`] and lives beside the format.
//!
//! # Every value has one home
//!
//! A domain's form holds what the domain states outside its lists, and each entry's form holds the
//! entry. Together they are exactly what [`crate::editable`] lists for the domain, no more and no
//! fewer — `every_value_has_one_home` checks that over every shipped scene. A form that left a
//! value out would be a panel that silently cannot express part of the file it is editing.
//!
//! # A menu's choices are the format's own
//!
//! A string the format reads as one of a set — a face, an axis, a schedule, how a face is held —
//! gets a menu, and the set is asked of the format rather than written here. The value is replaced
//! with a word no enum has and the scene parsed: serde's refusal names every variant it expected.
//! Each is then tried in turn, and only those the scene still parses with are offered, which is
//! what keeps a tag like a room's `"as": "mode"` from offering `pulse`, a word that would make the
//! keys beside it unknown. A key whose every string parses is free text, and `material` keeps the
//! catalogue the inspector always offered.

use crate::edit::{
    self, collect, escape, value_span, Choices, Editable, Value, IDENTITY, STRUCTURAL,
};
use pantometry_world::Scene;

/// The scene's own structure: what it states, before anything has run.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Model {
    /// The scene-level maps whose values are objects — `materials`, `composites`, `poses` — in
    /// file order, each with its keys.
    pub groups: Vec<Listing>,
    /// Every domain the scene states, in file order, whether or not it has an extent.
    pub domains: Vec<ModelDomain>,
}

/// A named set of entries: a scene-level map, or a list inside a domain.
#[derive(Clone, Debug, PartialEq)]
pub struct Listing {
    /// The key, as the file spells it.
    pub key: String,
    /// What each entry is called in the tree, in file order. A map's are its keys; a list's are
    /// the entry's `name` or `label`, else the strings it states, else the key and a number.
    pub entries: Vec<String>,
}

/// One domain the scene states, and the lists it holds.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelDomain {
    /// Its `name`.
    pub name: String,
    /// Its `kind`, as the file spells it — shown, never dispatched on.
    pub kind: String,
    /// Its lists of entries, in file order.
    pub lists: Vec<Listing>,
}

/// The model tree of `text`. Empty when the text is not JSON: there is nothing to point at, and
/// the source pane is where that gets fixed.
pub fn model(text: &str) -> Model {
    let Ok(root) = serde_json::from_str::<serde_json::Value>(text) else {
        return Model::default();
    };
    let Some(object) = root.as_object() else {
        return Model::default();
    };
    let mut groups = Vec::new();
    for key in in_file_order(text, "", object.keys()) {
        let member = &object[&key];
        if key != "domains" && is_group(member) {
            let base = format!("/{}", escape(&key));
            let names = member
                .as_object()
                .map(|m| in_file_order(text, &base, m.keys()))
                .unwrap_or_default();
            groups.push(Listing {
                key,
                entries: names,
            });
        }
    }
    let mut domains = Vec::new();
    let listed = root.get("domains").and_then(|d| d.as_array());
    for (i, d) in listed.into_iter().flatten().enumerate() {
        // A domain with no name has no path to be selected by, and the check already refuses it.
        let Some(name) = d.get("name").and_then(|n| n.as_str()) else {
            continue;
        };
        let kind = d
            .get("kind")
            .and_then(|k| k.as_str())
            .unwrap_or("")
            .to_string();
        let base = format!("/domains/{i}");
        let mut lists = Vec::new();
        if let Some(members) = d.as_object() {
            for key in in_file_order(text, &base, members.keys()) {
                let member = &members[&key];
                if !is_list(&kind, &key, member) {
                    continue;
                }
                let at = format!("{base}/{}", escape(&key));
                let entries = member
                    .as_array()
                    .into_iter()
                    .flatten()
                    .enumerate()
                    .map(|(j, x)| entry_name(text, &format!("{at}/{j}"), x, &key, j))
                    .collect();
                lists.push(Listing { key, entries });
            }
        }
        domains.push(ModelDomain {
            name: name.to_string(),
            kind,
            lists,
        });
    }
    Model { groups, domains }
}

/// Keys of the object at `base`, in the order the file writes them.
///
/// `serde_json` here backs an object with a `BTreeMap`, so its keys come out alphabetised — see
/// [`crate::edit`]. A tree and a form in the author's order read the way the file does: `kind` and
/// `name` first, a face before the coefficients that describe it.
fn in_file_order<'a>(
    text: &str,
    base: &str,
    keys: impl Iterator<Item = &'a String>,
) -> Vec<String> {
    let mut keyed: Vec<(usize, String)> = keys
        .map(|k| {
            let at =
                value_span(text, &format!("{base}/{}", escape(k))).map_or(usize::MAX, |r| r.start);
            (at, k.clone())
        })
        .collect();
    keyed.sort();
    keyed.into_iter().map(|(_, k)| k).collect()
}

/// A scene-level map that is a set of things — every value an object, and at least one.
fn is_group(member: &serde_json::Value) -> bool {
    member
        .as_object()
        .is_some_and(|m| !m.is_empty() && m.values().all(|v| v.is_object()))
}

/// A domain's member that is a list of entries: every element an object, or no elements and a
/// starting entry the format keeps for it — an empty `[]` is a list somebody emptied, and it can
/// be added to again.
fn is_list(kind: &str, key: &str, member: &serde_json::Value) -> bool {
    match member.as_array() {
        Some(items) if items.is_empty() => pantometry_world::templates::entry(kind, key).is_some(),
        Some(items) => items.iter().all(|x| x.is_object()),
        None => false,
    }
}

/// What an entry is called in the tree.
///
/// Its `name` or `label` when it has one, since that is what other entries call it. Otherwise the
/// strings it states, in file order — a face, a material, the two compartments a link joins —
/// which is what tells two cooled faces apart. An entry stating nothing but numbers is its list's
/// key and its place in it.
fn entry_name(text: &str, pointer: &str, value: &serde_json::Value, key: &str, j: usize) -> String {
    if let Some(id) = IDENTITY
        .iter()
        .find_map(|k| value.get(*k).and_then(|v| v.as_str()))
    {
        return id.to_string();
    }
    let mut said = Vec::new();
    if let Some(members) = value.as_object() {
        for k in in_file_order(text, pointer, members.keys()) {
            match &members[&k] {
                serde_json::Value::String(s) => said.push(s.clone()),
                serde_json::Value::Array(items)
                    if !items.is_empty() && items.iter().all(|x| x.is_string()) =>
                {
                    let words: Vec<&str> = items.iter().filter_map(|x| x.as_str()).collect();
                    said.push(words.join(", "));
                }
                _ => {}
            }
        }
    }
    if said.is_empty() {
        format!("{key} {}", j + 1)
    } else {
        said.join(", ")
    }
}

/// What the settings pane shows for one row of the tree, and what it can do there.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Form {
    /// The values this row states, in file order, each with the widget it wants.
    pub rows: Vec<Row>,
    /// The domain this row is, or is part of.
    pub domain: Option<String>,
    /// The lists of [`Form::domain`] that **add** can start an entry in: for a domain row, each
    /// list its kind has a starting entry for and the domain does not hold yet; for a list row,
    /// that list.
    pub adds: Vec<String>,
    /// For an entry row, which list it is in and where, so it can be duplicated and removed.
    pub entry: Option<(String, usize)>,
}

/// One line of a form.
#[derive(Clone, Debug, PartialEq)]
pub enum Row {
    /// One value.
    One(Editable),
    /// The elements of one key's array — a cell count, a corner, a position — side by side.
    Many {
        /// The key, as the file spells it.
        label: String,
        /// The unit the key names, or empty.
        unit: &'static str,
        /// The elements, in order.
        items: Vec<Editable>,
    },
}

impl Form {
    /// Every value the form offers, in order, whichever row it is drawn in.
    pub fn values(&self) -> Vec<&Editable> {
        self.rows
            .iter()
            .flat_map(|r| match r {
                Row::One(e) => vec![e],
                Row::Many { items, .. } => items.iter().collect(),
            })
            .collect()
    }
}

/// The form for the tree row at `path`.
///
/// `path` is a [`crate::Node::path`]. A domain's results rows, `/run/<name>` and
/// `/readings/<name>`, open the domain's form, because after a run those are the rows a person is
/// looking at. A row that states nothing — a group, a reading — has an empty form.
pub fn form(text: &str, path: &str) -> Form {
    let Ok(root) = serde_json::from_str::<serde_json::Value>(text) else {
        return Form::default();
    };
    let choices = Choices::of(&root);
    let mut out = Vec::new();
    let mut form = Form::default();
    let parts: Vec<&str> = path.trim_start_matches('/').split('/').collect();

    match parts.as_slice() {
        [""] => {
            // The scene itself: what it says outside its domains and its groups, each of which
            // has rows of its own.
            if let Some(object) = root.as_object() {
                for (key, member) in object {
                    if key == "domains" || STRUCTURAL.contains(&key.as_str()) || is_group(member) {
                        continue;
                    }
                    collect(
                        member,
                        &format!("/{}", escape(key)),
                        key,
                        &choices,
                        &mut out,
                    );
                }
            }
        }
        // **The rest of the path is the name, slashes and all.** A material may be called
        // `FR-4/Cu` — the format refuses only an empty name — and a group's entries have no rows
        // under them, so nothing after the key can be anything but the name. Split on every `/`,
        // that row's form was empty and its values had no home anywhere.
        ["scene", key, rest @ ..] if !rest.is_empty() => {
            let name = rest.join("/");
            if let Some(member) = root.get(*key).and_then(|g| g.get(&name)) {
                let at = format!("/{}/{}", escape(key), escape(&name));
                collect(member, &at, "", &choices, &mut out);
            }
        }
        ["domains" | "run" | "readings", name] => {
            if let Ok((i, d, kind)) = domain(&root, name) {
                if let Some(members) = d.as_object() {
                    for (key, member) in members {
                        if STRUCTURAL.contains(&key.as_str()) || is_list(kind, key, member) {
                            continue;
                        }
                        let at = format!("/domains/{i}/{}", escape(key));
                        collect(member, &at, key, &choices, &mut out);
                    }
                }
                form.domain = Some(name.to_string());
                // Only the lists the domain does not hold yet. One it holds has a row of its own,
                // and that row's add is the one place to add to it: two buttons for one command
                // is two answers to one question, which is what `editor-core`'s README says a
                // toolbar beside a menu was.
                form.adds = pantometry_world::templates::ENTRIES
                    .iter()
                    .filter(|e| {
                        e.kind == kind && !d.get(e.key).is_some_and(|m| is_list(kind, e.key, m))
                    })
                    .map(|e| e.key.to_string())
                    .collect();
            }
        }
        ["domains", name, key] => {
            if let Ok((_, d, kind)) = domain(&root, name) {
                if d.get(*key).is_some_and(|m| is_list(kind, key, m)) {
                    form.domain = Some(name.to_string());
                    if pantometry_world::templates::entry(kind, key).is_some() {
                        form.adds = vec![key.to_string()];
                    }
                }
            }
        }
        ["domains", name, key, j] => {
            if let (Ok((i, d, kind)), Ok(j)) = (domain(&root, name), j.parse::<usize>()) {
                let list = d.get(*key).filter(|m| is_list(kind, key, m));
                if let Some(entry) = list.and_then(|l| l.get(j)) {
                    let at = format!("/domains/{i}/{}/{j}", escape(key));
                    collect(entry, &at, "", &choices, &mut out);
                    form.domain = Some(name.to_string());
                    form.entry = Some((key.to_string(), j));
                }
            }
        }
        _ => {}
    }

    // Whether the scene parses as it stands, which decides how a menu's choices are tested.
    let parses = serde_json::from_str::<Scene>(text).is_ok();
    for field in &mut out {
        if let Value::Text { now, choices } = &mut field.value {
            if choices.is_empty() {
                *choices = format_choices(text, &field.pointer, now, parses);
            }
        }
    }
    out.sort_by_key(|f| value_span(text, &f.pointer).map_or(usize::MAX, |r| r.start));
    form.rows = rows(out);
    form
}

/// Domain `name`: where it is in `domains`, the object, and its kind.
fn domain<'a>(
    root: &'a serde_json::Value,
    name: &str,
) -> Result<(usize, &'a serde_json::Value, &'a str), ()> {
    let domains = root.get("domains").and_then(|d| d.as_array()).ok_or(())?;
    let i = domains
        .iter()
        .position(|d| d.get("name").and_then(|n| n.as_str()) == Some(name))
        .ok_or(())?;
    let kind = domains[i]
        .get("kind")
        .and_then(|k| k.as_str())
        .unwrap_or("");
    Ok((i, &domains[i], kind))
}

/// A word no enum in the format has, written where a string is to ask what the format expects.
const PROBE: &str = "pantometry-editor-probe";

/// The words the format accepts for the string at `pointer`, or empty when it takes any string.
///
/// See the module documentation for how: the refusal of a word nobody uses names the set, and each
/// member of the set is kept only if the scene still parses with it. When the scene does not parse
/// as it stands that last step cannot tell anything, so the set is offered whole.
///
/// A set of one — a tag whose every other variant needs other keys — is no choice, and is left as
/// free text, which is what it was before.
fn format_choices(text: &str, pointer: &str, now: &str, parses: bool) -> Vec<String> {
    let Ok(probed) = edit::set_text(text, pointer, PROBE) else {
        return Vec::new();
    };
    let Err(refusal) = serde_json::from_str::<Scene>(&probed) else {
        return Vec::new();
    };
    let said = refusal.to_string();
    // The refusal has to be about the probe. A scene that does not parse for some other reason
    // can refuse at an earlier byte, and its message is not a list of anything here.
    let marker = format!("unknown variant `{PROBE}`");
    let Some(at) = said.find(&marker) else {
        return Vec::new();
    };
    let tail = &said[at + marker.len()..];
    let Some(expected) = tail.find("expected") else {
        return Vec::new();
    };
    let named: Vec<String> = tail[expected..]
        .split('`')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect();
    let kept: Vec<String> = named
        .into_iter()
        .filter(|word| {
            word == now
                || !parses
                || edit::set_text(text, pointer, word)
                    .is_ok_and(|t| serde_json::from_str::<Scene>(&t).is_ok())
        })
        .collect();
    if kept.len() < 2 {
        Vec::new()
    } else {
        kept
    }
}

/// Fields in file order, with the elements of one array gathered into one row.
fn rows(fields: Vec<Editable>) -> Vec<Row> {
    let mut out: Vec<Row> = Vec::new();
    for field in fields {
        if let Some((base, index)) = element_of(&field.label) {
            if let Some(Row::Many { label, items, .. }) = out.last_mut() {
                if *label == base && items.len() == index {
                    items.push(field);
                    continue;
                }
            }
            if index == 0 {
                out.push(Row::Many {
                    label: base,
                    unit: field.unit,
                    items: vec![field],
                });
                continue;
            }
        }
        out.push(Row::One(field));
    }
    out
}

/// `cells[2]` is element 2 of `cells`. `regions[0].material` is not an element: it is a value
/// inside one, and is drawn on a line of its own.
fn element_of(label: &str) -> Option<(String, usize)> {
    let open = label.strip_suffix(']')?.rfind('[')?;
    let index = label[open + 1..label.len() - 1].parse().ok()?;
    Some((label[..open].to_string(), index))
}
