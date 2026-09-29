//! Layered world memory for embodied agents.
//!
//! This is external state, not a Transformer KV cache. Static geometry, persistent
//! movable objects, and live dynamic tracks have different update rules and lifetimes.
//! A model or planner can query the relevant local slice without re-encoding the whole
//! environment on every step.

use crate::geo::Enu;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PersistenceClass {
    /// Structural geometry expected to remain fixed: walls, doorsills, columns.
    Structural,
    /// Usually stable but physically movable: tables, chairs, cabinets.
    PersistentMovable,
    /// Short-lived or continuously moving: people, animals, carried objects.
    Dynamic,
}

#[derive(Clone, Debug, PartialEq)]
pub struct WorldObject {
    pub id: String,
    pub label: Option<String>,
    pub class: PersistenceClass,
    pub pos: Enu,
    pub vel: Enu,
    pub confidence: f64,
    /// Age of the latest evidence supporting this state.
    pub evidence_age_s: f64,
    /// Monotonic revision number. Consumers can request only objects changed since a revision.
    pub revision: u64,
}

#[derive(Clone, Debug, Default)]
pub struct WorldMemory {
    structural: HashMap<String, WorldObject>,
    persistent: HashMap<String, WorldObject>,
    dynamic: HashMap<String, WorldObject>,
    revision: u64,
}

impl WorldMemory {
    pub fn revision(&self) -> u64 {
        self.revision
    }

    pub fn upsert(&mut self, mut object: WorldObject) {
        self.revision = self.revision.saturating_add(1);
        object.revision = self.revision;

        // An id must exist in exactly one layer. Classification changes are explicit.
        self.structural.remove(&object.id);
        self.persistent.remove(&object.id);
        self.dynamic.remove(&object.id);

        match object.class {
            PersistenceClass::Structural => {
                self.structural.insert(object.id.clone(), object);
            }
            PersistenceClass::PersistentMovable => {
                self.persistent.insert(object.id.clone(), object);
            }
            PersistenceClass::Dynamic => {
                self.dynamic.insert(object.id.clone(), object);
            }
        }
    }

    pub fn get(&self, id: &str) -> Option<&WorldObject> {
        self.dynamic
            .get(id)
            .or_else(|| self.persistent.get(id))
            .or_else(|| self.structural.get(id))
    }

    /// Objects changed after a consumer's last known revision.
    pub fn changed_since(&self, revision: u64) -> Vec<&WorldObject> {
        self.structural
            .values()
            .chain(self.persistent.values())
            .chain(self.dynamic.values())
            .filter(|o| o.revision > revision)
            .collect()
    }

    /// Drop stale dynamic state without touching the static map or persistent furniture.
    pub fn prune_dynamic(&mut self, max_age_s: f64) {
        self.dynamic
            .retain(|_, o| o.evidence_age_s.is_finite() && o.evidence_age_s <= max_age_s);
    }

    pub fn structural(&self) -> impl Iterator<Item = &WorldObject> {
        self.structural.values()
    }

    pub fn persistent(&self) -> impl Iterator<Item = &WorldObject> {
        self.persistent.values()
    }

    pub fn dynamic(&self) -> impl Iterator<Item = &WorldObject> {
        self.dynamic.values()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obj(id: &str, class: PersistenceClass) -> WorldObject {
        WorldObject {
            id: id.into(),
            label: None,
            class,
            pos: Enu { e: 0.0, n: 0.0, u: 0.0 },
            vel: Enu { e: 0.0, n: 0.0, u: 0.0 },
            confidence: 1.0,
            evidence_age_s: 0.0,
            revision: 0,
        }
    }

    #[test]
    fn static_geometry_survives_dynamic_pruning() {
        let mut memory = WorldMemory::default();
        memory.upsert(obj("wall-1", PersistenceClass::Structural));

        let mut person = obj("person-1", PersistenceClass::Dynamic);
        person.evidence_age_s = 5.0;
        memory.upsert(person);

        memory.prune_dynamic(1.0);

        assert!(memory.get("wall-1").is_some());
        assert!(memory.get("person-1").is_none());
    }

    #[test]
    fn changed_since_returns_only_new_revisions() {
        let mut memory = WorldMemory::default();
        memory.upsert(obj("wall-1", PersistenceClass::Structural));
        let rev = memory.revision();
        memory.upsert(obj("glass-1", PersistenceClass::Dynamic));

        let changed = memory.changed_since(rev);
        assert_eq!(changed.len(), 1);
        assert_eq!(changed[0].id, "glass-1");
    }

    #[test]
    fn reclassification_moves_object_between_layers() {
        let mut memory = WorldMemory::default();
        memory.upsert(obj("chair-1", PersistenceClass::Dynamic));
        memory.upsert(obj("chair-1", PersistenceClass::PersistentMovable));

        assert_eq!(memory.dynamic().count(), 0);
        assert_eq!(memory.persistent().count(), 1);
    }
}
